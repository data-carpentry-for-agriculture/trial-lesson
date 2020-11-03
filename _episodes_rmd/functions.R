### These are functions found in the Data Carpentry for Agriculture Workshop

fitted_line <- function(data, x, y, x_label, y_label){
  ggplot() +
  geom_smooth(data = data, method = "gam", aes(y = y, x = x), size = 0.5, se=FALSE) +
  ylab(y_label) +
  xlab(x_label) +
  theme_grey(base_size = 12)
}

deposit_on_grid <- function(grid, data, col, fn = median){
  grid_sp <- as(grid, "Spatial")
  merge <- sp::over(as(grid, "Spatial"), as(data[, col], "Spatial"), fn = fn)
  grid_sp@data <- cbind(merge, grid_sp@data)
  subplots_data <- st_as_sf(grid_sp)
  return(subplots_data)
}

utm_zone <- function(long){
  utm <- (floor((long + 180)/6) %% 60) + 1
  return(utm)
}

st_transform_utm <- function(sfobject){
  crs <- st_crs(sfobject)
  if (str_detect(crs$wkt, "longitude") != TRUE){
    print("Not in lat/long. Returning original object.")
    return(sfobject)
  }
  else {
    utmzone <- utm_zone(mean(st_bbox(sfobject)[c(1,3)]))
    projutm <- as.numeric(paste0("326", utmzone))
    newobj <- st_transform(sfobject, projutm)
    return(newobj)
  }
}

# `long2UTM()` is a function written to take the argument `long` and output the
# result of the equation. `floor()` returns the largest integer that is not
# greater than the input.  The `st_bbox()` function returns the bounding box of
# the dataset, i.e. the four corners defining a rectangle that would contain all
# of the polygons.  By taking the mean of the first and third items returned by
# `st_bbox()`, we get the longitude of the point directly in the center of the box.
#
# The code below pastes together the full
# ESPG code for any `utmzone` we calculate with `long2UTM`. `paste0()` pastes
# together the two arguments `"326"` and `utmzone` as string. But when we transform
# `trial` into UTM, we only need the ESPG number. So we convert that to numeric with
# `as.numeric()`, giving us a final ESPG of `r paste0("326", utmzone)`.

cent_long <- function(sfobject){
  mean(st_bbox(sfobject)[c(1,3)])
}

cent_lat <- function(sfobject){
  mean(st_bbox(sfobject)[c(2,4)])
}

as.Date.daymetr <- function(daymetdata){
  date <- as.Date(daymetdata$yday, origin = paste0(daymetdata$year-1, "-12-31"))
  return(date)
}

sumprec_by_monthyear <- function(data){
  data %>% dplyr::group_by(month, year) %>% dplyr::summarise(prec_month = sum(prec))
}

avgprec_by_month <- function(data){
  data %>% dplyr::group_by(month) %>% dplyr::summarise(prec_avg = mean(prec_month))
}

# conversions
c_to_f <- function(varc){
  conv_unit(varc, "C", "F")
}

mm_to_in <- function(varc){
  conv_unit(varc, "mm", "inch")
}

cm_to_in <- function(varc){
  conv_unit(varc, "cm", "inch")
}

ha_to_ac <- function(varc){
  conv_unit(varc, "hectare", "acres")
}

kg_to_lb <- function(varc){
  conv_unit(varc, "kg", "lb")
}

m_to_ft <- function(varc){
  conv_unit(varc, "m", "ft")
}

map_poly <- function(sfobject, variable, name, colors="default"){
  if (colors != "default"){
     tm_shape(sfobject) + tm_polygons(variable, title = name) +
       tm_layout(legend.outside = TRUE, frame = FALSE, palette=colors)
  } else {
     tm_shape(sfobject) + tm_polygons(variable, title = name) +
       tm_layout(legend.outside = TRUE, frame = FALSE)
  }
}

# in case we need to see any plots better we can mess with the colors: https://geocompr.robinlovelace.net/adv-map.html
map_points <- function(sfobject, variable, name, colors="default"){
  if (colors != "default"){
    tm_shape(sfobject) + tm_dots(variable, title = name, palette=colors) +
    tm_layout(legend.outside = TRUE, frame = FALSE)
  } else {
    tm_shape(sfobject) + tm_dots(variable, title = name) +
      tm_layout(legend.outside = TRUE, frame = FALSE)
      }
}

wt_mean <- function(property, weights)
{
  # compute thickness weighted mean, but only when we have enough data
  # in that case return NA

  # save indices of data that is there
  property.that.is.na <- which( is.na(property) )
  property.that.is.not.na <- which( !is.na(property) )

  if( length(property) - length(property.that.is.na) >= 1)
    prop.aggregated <- sum(weights[property.that.is.not.na] * property[property.that.is.not.na], na.rm=TRUE) / sum(weights[property.that.is.not.na], na.rm=TRUE)
  else
    prop.aggregated <- NA

  return(prop.aggregated)
}

profile_total <- function(property, thickness)
{
  # compute profile total
  # in that case return NA

  # save indices of data that is there
  property.that.is.na <- which( is.na(property) )
  property.that.is.not.na <- which( !is.na(property) )

  if( length(property) - length(property.that.is.na) >= 1)
    prop.aggregated <- sum(thickness[property.that.is.not.na] * property[property.that.is.not.na], na.rm=TRUE)
  else
    prop.aggregated <- NA

  return(prop.aggregated)
}

# define a function to perfom hz-thickness weighted aggregtion
component_level_aggregation <- function(i)
{
  # horizon thickness is our weighting vector
  hz_thick <- i$hzdepb.r - i$hzdept.r

  # compute wt.mean aggregate values
  clay <- wt_mean(i$claytotal.r, hz_thick)
  silt <- wt_mean(i$silttotal.r, hz_thick)
  sand <- wt_mean(i$sandtotal.r, hz_thick)
  # compute profile sum values
  water_storage <- profile_total(i$awc.r, hz_thick)

  # make a new dataframe out of the aggregate values
  d <- data.frame(cokey=unique(i$cokey), clay=clay, silt=silt, sand=sand, water_storage=water_storage)

  return(d)
}

mapunit_level_aggregation <- function(i)
{
  # component percentage is our weighting vector
  comppct <- i$comppct.r

  # wt. mean by component percent
  clay <- wt_mean(i$clay, comppct)
  silt <- wt_mean(i$silt, comppct)
  sand <- wt_mean(i$sand, comppct)
  water_storage <- wt_mean(i$water_storage, comppct)

  # make a new dataframe out of the aggregate values
  d <- data.frame(mukey=unique(i$mukey), clay=clay, silt=silt, sand=sand, water_storage=water_storage)

  return(d)
}

###############################################
# using these functions on my data
c_s_s_soil <- function(ssurgo){

  chorizon<-ssurgo$tabular$chorizon

  component<-ssurgo$tabular$component

  # aggregate horizon data to the component level
  chorizon.agg <- ddply(chorizon, .(cokey), .fun=component_level_aggregation, .progress='text')

  # join up the aggregate chorizon data to the component table
  comp.merged <- merge(component, chorizon.agg, by='cokey')

  # aggregate component data to the map unit level
  component.agg <- ddply(comp.merged, .(mukey), .fun=mapunit_level_aggregation, .progress='text')

  # save data to the ssurgo shapefile
  spatial <- as(ssurgo$spatial, "sf")
  spatial <- dplyr::rename(spatial, mukey = MUKEY)
  spatial <- merge.data.frame(spatial, component.agg, by="mukey")
  spatial <- spatial[,-5]
  spatial <- dplyr::distinct(spatial)
  #spatial <- spatial[,6:9]
  return(spatial)
}

clean_sd <- function(data, var, sd_no){
  data <- subset(data, var >= mean(var, na.rm = TRUE) - sd_no*sd(var, na.rm = TRUE) &  var <= mean(var, na.rm = TRUE) + sd_no*sd(var, na.rm = TRUE))
  return(data)
}

make_abline <- function(LongA,LongB,LatA,LatB,projutm){
  ab_string <- rbind(c(LongA, LatA),c(LongB, LatB)) %>%
    st_linestring()
  ab_line <-st_sf(
    id='ab_line',
    st_sfc(ab_string)
  )
  st_crs(ab_line) <- '+proj=longlat +datum=WGS84 +no_defs'
  ab_line <- st_transform(ab_line,projutm)
  return(ab_line)
}

make_grids <- function(trialarea, ab_line, long_in, short_in, length_ft, width_ft, set_seed=FALSE){

  # for the purposes of the workshop
  if (set_seed){
    set.seed(42) # 42, obvs :D
  }
  bbox_field <- st_bbox(trialarea)
  # +c(100,100,50,50)
  xmin <- bbox_field[1]
  ymin <- bbox_field[2]
  xmax <- bbox_field[3]
  ymax <- bbox_field[4]

  #--- identify the starting point ---#
  starting_point <- c(xmax,ymax)

  plot_length_meter <- conv_unit(length_ft, 'ft', 'm')
  plot_width_meter <- conv_unit(width_ft, 'ft', 'm')

  #===================================
  # Define the long and short vectors
  #===================================
  ab_1 <- st_geometry(ab_line)[[1]][1,]
  ab_2 <- st_geometry(ab_line)[[1]][2,]

  #--------------------------
  # find the origin, end point, and the rotation matrix
  #--------------------------
  if (long_in=='SN'){ # if the plot is long in SN direction
    #--- find the origin and end points ---#
    if (ab_1[2]>=ab_2[2]){
      origin <- ab_2
      end_point <- ab_1
    } else{
      origin <- ab_1
      end_point <- ab_2
    }
    #--- find rotation vector ---#
    if (short_in=='WE'){
      rotate_mat <- matrix(c(cos(-3.000*pi/2),sin(-3.000*pi/2),-sin(-3.000*pi/2),cos(-3.000*pi/2)),nrow=2)
    } else{
      rotate_mat <- matrix(c(cos(3.000*pi/2),sin(3.000*pi/2),-sin(3.000*pi/2),cos(3.000*pi/2)),nrow=2)
    }
  } else if (long_in=='NS') {
    #--- find the origin and end points ---#
    if (ab_1[2]>=ab_2[2]){
      origin <- ab_1
      end_point <- ab_2
    } else{
      origin <- ab_2
      end_point <- ab_1
    }
    #--- find rotation vector ---#
    if (short_in=='WE'){
      rotate_mat <- matrix(c(cos(3.000*pi/2),sin(3.000*pi/2),-sin(3.000*pi/2),cos(3.000*pi/2)),nrow=2)
    } else{
      rotate_mat <- matrix(c(cos(-3.000*pi/2),sin(-3.000*pi/2),-sin(-3.000*pi/2),cos(-3.000*pi/2)),nrow=2)
    }
  } else if (long_in=='WE') {
    #--- find the origin and end points ---#
    if (ab_1[1]>=ab_2[1]){
      origin <- ab_2
      end_point <- ab_1
    } else{
      origin <- ab_1
      end_point <- ab_2
    }
    #--- find rotation vector ---#
    if (short_in=='SN'){
      rotate_mat <- matrix(c(cos(-3.000*pi/2),sin(-3.000*pi/2),-sin(-3.000*pi/2),cos(-3.000*pi/2)),nrow=2)
    } else if (short_in=='NS'){
      rotate_mat <- matrix(c(cos(3.000*pi/2),sin(3.000*pi/2),-sin(3.000*pi/2),cos(3.000*pi/2)),nrow=2)
    }
  } else if (long_in=='EW'){
    #--- find the origin and end points ---#
    if (ab_1[1]>=ab_2[1]){
      origin <- ab_1
      end_point <- ab_2
    } else{
      origin <- ab_2
      end_point <- ab_1
    }

    #--- find rotation vector ---#
    if (short_in=='SN'){
      rotate_mat <- matrix(c(cos(3.000*pi/2),sin(3.000*pi/2),-sin(3.000*pi/2),cos(3.000*pi/2)),nrow=2)
    } else if (short_in=='NS'){
      rotate_mat <- matrix(c(cos(-3.000*pi/2),sin(-3.000*pi/2),-sin(-3.000*pi/2),cos(-3.000*pi/2)),nrow=2)
    }
  }

  #--------------------------
  # Find the long and short vectors
  #--------------------------
  #--- long vector ---#
  long_vec <- end_point - origin

  #--- short vector ---#
  short_vec <- rotate_mat %*% long_vec

  #--------------------------
  # normalize the vectors
  #--------------------------
  vector_len <- sqrt(long_vec[1]^2+long_vec[2]^2)
  long_norm <- long_vec/vector_len*plot_length_meter
  short_norm <- (short_vec/vector_len*plot_width_meter) %>% as.vector()

  #===================================
  # Create grids
  #===================================

  bbox_field['ymax']-bbox_field['ymin']



  #--- how many rows and columns ---#
  if (long_in %in% c('SN','NS')){
    num_rows <- ceiling((bbox_field['ymax']-bbox_field['ymin'])/plot_length_meter)
    num_cols <- ceiling((bbox_field['xmax']-bbox_field['xmin'])/plot_width_meter) + 1 #want an extra to cover borders.
  } else if (long_in %in% c('WE','EW')){
    num_rows <- ceiling((bbox_field['ymax']-bbox_field['ymin'])/plot_width_meter)
    num_cols <- ceiling((bbox_field['xmax']-bbox_field['xmin'])/plot_length_meter) + 1 #want an extra to cover borders.
  }

  #--------------------------
  # Create grids
  #--------------------------

  all_polygons_ls <- list()
  if (long_in %in% c('SN','NS')){ # if the applicator moves NS or SN
    if (short_in %in% c('EW')){
      for (i in 1:num_cols){
        #i=1
        if(i==1){
          col_start <- starting_point
        } else{
          col_start <- st_geometry(all_polygons_ls[[i-1]])[[1]][[1]][4,]
        }

        col_polygons_ls <- list()
        for (j in 1:num_rows){
          #j=1
          if(j==1){
            point_1 <- col_start
          } else{
            point_1 <- col_polygons_ls[[j-1]][[1]][2,]
          }

          point_2 <- point_1+long_norm
          point_3 <- point_2-short_norm
          point_4 <- point_3-long_norm

          p_temp <- rbind(point_1,point_2,point_3,point_4,point_1) %>%
            list() %>%
            st_polygon()
          col_polygons_ls[[j]] <- p_temp
        }


        all_polygons_ls[[i]] <-
          st_sf(
            plotid=(1+(i-1)*num_rows):(i*num_rows),
            GRIDY = (1:num_rows),
            GRIDX = i,
            st_sfc(col_polygons_ls)
          )
      }
    }else{#needs to complete
    }

  }else if(long_in %in% c('WE','EW')){ # if the applicator moves WE or EW
    if (short_in %in% c('SN')){
      for (i in 1:num_rows){
        # i=2
        if(i==1){
          col_start <- starting_point
        } else{
          col_start <- st_geometry(all_polygons_ls[[i-1]])[[1]][[1]][2,]
        }

        col_polygons_ls <- list()
        for (j in 1:num_cols){
          #j = 2
          if(j==1){
            point_1 <- col_start
          } else{
            point_1 <- col_polygons_ls[[j-1]][[1]][4,]
          }
          point_2 <- point_1+short_norm
          point_3 <- point_2+long_norm
          point_4 <- point_3-short_norm

          # point_2 <- point_1+long_norm
          # point_3 <- point_2+short_norm
          # point_4 <- point_3-long_norm

          p_temp <- rbind(point_1,point_2,point_3,point_4,point_1) %>%
            list() %>%
            st_polygon()
          col_polygons_ls[[j]] <- p_temp
        }

        all_polygons_ls[[i]] <-
          st_sf(
            plotid=(1+(i-1)*num_cols):(i*num_cols),
            GRIDY = i,
            GRIDX = (1:num_cols),
            st_sfc(col_polygons_ls)
          )
      }
    }else if(short_in %in% c('NS')){
      for (i in 1:num_rows){
        # i=2
        if(i==1){
          col_start <- starting_point
        } else{
          col_start <- st_geometry(all_polygons_ls[[i-1]])[[1]][[1]][1,]
        }

        col_polygons_ls <- list()
        for (j in 1:num_cols){
          #j = 2
          if(j==1){
            point_2 <- col_start
          } else{
            point_2 <- col_polygons_ls[[j-1]][[1]][3,]
          }
          point_1 <- point_2+short_norm
          point_3 <- point_2+long_norm
          point_4 <- point_3+short_norm

          # point_2 <- point_1+long_norm
          # point_3 <- point_2+short_norm
          # point_4 <- point_3-long_norm

          p_temp <- rbind(point_1,point_2,point_3,point_4,point_1) %>%
            list() %>%
            st_polygon()
          col_polygons_ls[[j]] <- p_temp
        }

        all_polygons_ls[[i]] <-
          st_sf(
            plotid=(1+(i-1)*num_cols):(i*num_cols),
            GRIDY = i,
            GRIDX = (1:num_cols),
            st_sfc(col_polygons_ls)
          )
      }
    }
  }
  #--- combine all the grids ---#
  all_grids <- do.call(rbind,all_polygons_ls)
  all_grids <- dplyr::rename(all_grids, geom = st_sfc.col_polygons_ls.)
  return(all_grids)
}

make_subplots <- function(boundary.utm,ab_line,long_in,short_in,starting_point){

  design_grids_utm <- make_grids(boundary.utm,ab_line,long_in,short_in,starting_point)
  st_crs(design_grids_utm) <-st_crs(boundary.utm) #bothfields file came with crs as utm

  CentX <- lapply(st_geometry(st_centroid(design_grids_utm)),function(x) x[1]) %>% unlist()
  CentY <- lapply(st_geometry(st_centroid(design_grids_utm)),function(x) x[2]) %>% unlist()
  design_grids_utm <- design_grids_utm %>%
    mutate(cent_x=CentX,cent_y=CentY)

  trial_grid1_utm <- st_intersection(boundary.utm, design_grids_utm)
  trial_grid1_utm$Name <- NULL
  rownames(trial_grid1_utm) <- 1:nrow(trial_grid1_utm)

  max_area <- as.numeric(max(st_area(trial_grid1_utm)))

  trial_grid <- trial_grid1_utm %>%
    mutate(area=as.numeric(st_area(.))) %>%
    mutate(drop=ifelse(area<(max_area-100), 1, 0) )

  return(trial_grid)
}

treat_assign <- function(trialarea, trial_grid, head_buffer_ft, seed_treat_rates, nitrogen_treat_rates, seed_quo, nitrogen_quo, set_seed=FALSE){
  if (set_seed){
    set.seed(42) # obvs :D
  }
  head_buffer_m <- conv_unit(head_buffer_ft, 'ft', 'm')
  infield <- st_buffer(trialarea, -head_buffer_m)
  outfield <- st_difference(trial_grid, infield)

  intrial <- st_intersection(trial_grid, infield)

  intrial$dummy <- 1
  outfield$dummy <- 0

  trial_grid <- rbind(intrial, outfield)

  max_area <- as.numeric(mean(st_area(trial_grid)))

  trial_grid <- trial_grid %>%
    mutate(area = as.numeric(st_area(.))) %>%
    mutate(small = ifelse(area < (max_area*0.9), 1, 0)) %>%
    mutate(drop = ifelse(dummy == 0 | small == 1, 1, 0))  # drop polygons in headlands from dummy
    # drop polygons that are too small

  tm_shape(trial_grid) + tm_polygons("drop")

  num_treats <- length(seed_treat_rates)*length(nitrogen_treat_rates)
  num_plots <- nrow(trial_grid)
  num_rep <- floor(num_plots/num_treats)

  treat_ls <- list()
  for (i in 1:num_rep){
    treat_ls[[i]] <- sample(1:num_treats,num_treats,replace=FALSE)
  }

  remainder_ls <- sample(1:num_treats,num_plots%%num_treats,replace=FALSE)
  treat_list <- c(unlist(treat_ls),remainder_ls)

  trial_grid$GRIDID <- trial_grid$plotid


  grid_list <- trial_grid %>%
    .$GRIDID

  grid_to_treat <- data.table(GRIDID=grid_list,treat_type=treat_list)

  trial_grid3 <-left_join(trial_grid,grid_to_treat,by='GRIDID') %>%
    mutate(treat_type=ifelse(drop == 1,num_treats+1, treat_type)) %>%
    dplyr::select(-area,-drop)

  exception <- data.table(
    SEEDRATE=c(seed_quo),
    NRATE=c(nitrogen_quo),
    treat_type=c(num_treats+1)
  )

  pair_ls <- expand.grid(seed_treat_rates, nitrogen_treat_rates) %>%
    data.table() %>%
    setnames(names(.),c('SEEDRATE', 'NRATE')) %>%
    .[,treat_type:=1:num_treats] %>%
    rbind(exception) %>%
    dplyr::select(NRATE,SEEDRATE,treat_type)

  whole <- left_join(trial_grid3, pair_ls, by = 'treat_type')
  class(whole)

  tm_shape(whole) + tm_polygons("NRATE")

  whole_td <- broom::tidy(as(whole, "Spatial"))

  whole_td$new_id <- paste('ID',whole_td$id,sep='')
  whole_td$id <- whole_td$new_id
  whole_td$new_id <- NULL

  temp_whole <- whole %>%
    mutate(id=paste('ID',1:nrow(whole),sep='')) %>%
    dplyr::select(id,treat_type,NRATE,SEEDRATE)

  tm_shape(temp_whole) + tm_polygons("NRATE")

  return(temp_whole)
  if (set_seed){ #reset
    set.seed(NULL)
  }
}

st_over <- function(x, y) {
  sapply(sf::st_intersects(x, y), function(z)
    if (length(z) == 0)
      NA_integer_
    else
      z[1])
}

download_ssurgo <- function(name_of_field, boundary_sp_in, redo=FALSE){
  tryCatch(
    expr = {
      ssurgo_ <- get_ssurgo(boundary_sp_in, name_of_field, force.redo=redo)
      message("Successfully downloaded SSURGO.")
      return(ssurgo_)
    },
    error = function(e){
      message('Caught an error!  Let\'s try again after waiting 5 seconds...')
      print(e)
      Sys.sleep(5.0)
      download_ssurgo(name_of_field, boundary_sp_in)
    },
    warning = function(w){
      message('Caught a warning!')
      print(w)
    }
  )
}

clean_buffer <- function(buffer_object, buffer_ft, data){
  buffer_m <- conv_unit(buffer_ft, "ft", "m")
  buffer <- st_buffer(buffer_object, -buffer_m) # plots are 24 m wide and 2 yield passes
  ov <- st_over(data, st_geometry(buffer))
  data$out <- is.na(ov) # demarcate the yield values removed
  clean <- subset(data, out == FALSE)
  return(clean)
}


# JPN: for simulating yields, asapplied & asplanted based on trial design
# NOTE: known issues with slowness of rbind: https://stackoverflow.com/questions/14693956/how-can-i-prevent-rbind-from-geting-really-slow-as-dataframe-grows-larger

simulate_trial <- function(whole_plot, yield, asapplied, asplanted, useGLM=TRUE){
 # "yield" is an input, "planting" and "nitrogen" since those have changed

 # replace trial data with whole plot
 #trial <- whole_plot
 #asapplied <- nitrogen
 #asplanted <- planting

 # from fit
 if (useGLM){
  coefs = read.csv('https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/data/coefs_fit_glm.csv')
 } else {
  coefs = read.csv('https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/data/coefs_fit.csv')
 }
 # transform if needed
 if (st_crs(yield) != st_crs(whole_plot)){
   yieldutm = st_transform_utm(yield)
 }
 if (st_crs(asplanted) != st_crs(whole_plot)){
   asplanted = st_transform_utm(asplanted)
 }
 if (st_crs(asapplied) != st_crs(whole_plot)){
   asapplied = st_transform_utm(asapplied)
 }

 # also, add in random bigs
 randomBigProb = 0.005 # will pull random big, looks like this happens ~0.003 of the time in original trial data
 maxBig = 1200
 minBig = 400

 randomBigProbApp = 0.005 # random as applied
 maxBigApp = 50
 minBigApp = 1300


 # another param, how often to print
 nPrint = 50


 # loop through each geometry
 flag = 0 # flag to turn off one
 flagapp = 0
 flaggplant = 0
 print("This might take a little while... now is a great time for a coffee :)")
 for (i in 1:length(whole_plot$geom)){
  myOut2 = 0
  asappliedOut2 = 0
  asplantedOut2 = 0
  samps = 0
 #for (i in 1:3){ # test
   if (i%%nPrint==0){
     print(paste0('On ', i, ' of ', length(whole_plot$geom), ' geometries'))
   }
   yield_int <- st_intersection(yieldutm, whole_plot$geom[i])
   asapplied_int <- st_intersection(asapplied, whole_plot$geom[i])
   asplanted_int <- st_intersection(asplanted, whole_plot$geom[i])
   if (flagapp == 0){
       asappliedOut = asapplied_int
       if (nrow(asapplied_int)>0){
	 asappliedOut$Rate_Appli = whole_plot$NRATE[i]
       }
       flagapp = 1
   } else {
     if (nrow(asapplied_int)>0){
       asappliedOut2 = asapplied_int
       asappliedOut2$Rate_Appli = whole_plot$NRATE[i]
       # random wrong applications
       samps = runif(length(asappliedOut2$Rate_Appli))
       asappliedOut2$Rate_Appli[samps <= randomBigProbApp] = samps[samps <= randomBigProbApp]/randomBigProbApp*(maxBigApp-minBigApp) + minBigApp
       asappliedOut = rbind(asappliedOut, asappliedOut2)
     }
   }
   if (flaggplant == 0){
     asplantedOut = asplanted_int
     if (nrow(asplanted_int)>0){
       asplantedOut$Rt_Apd_Ct_ = whole_plot$SEEDRATE[i]
     }
     flaggplant = 1
   } else {
     if (nrow(asplanted_int)>0){
       asplantedOut2 = asplanted_int
       asplantedOut2$Rt_Apd_Ct_ = whole_plot$SEEDRATE[i]
       asplantedOut = rbind(asplantedOut, asplantedOut2)
     }
   }


   if (length(row(yield_int)) > 0){ # have entries, update
     # grab random index of row for coefficients of fit
     if (nrow(asplanted_int)>0){
       ele = mean(asplanted_int$Elevation_)
     } else {
       ele = mean(asplanted$Elevation_)
     }
     mycoefs = coefs[sample(nrow(coefs), nrow(yield_int)), ]
     yieldsMod = mycoefs[,'X.Intercept.'] + mycoefs[, 'Rate_Appli']*whole_plot$NRATE[i] +
       mycoefs[, 'Rt_Apd_Ct_']*whole_plot$SEEDRATE[i] + mycoefs[, 'Elevation_']*ele
     ## add in big stuff randomly
     samps = runif(length(yieldsMod))
     yieldsMod[samps <= randomBigProb] = samps[samps <= randomBigProb]/randomBigProb*(maxBig-minBig) + minBig
     if (flag == 0){
       myOut = yield_int
       myOut$Yld_Vol_Dr = yieldsMod
       flag = 1
     } else {
       myOut2 = yield_int
       myOut2$Yld_Vol_Dr = yieldsMod
       myOut = rbind(myOut, myOut2)
     }
   } else { # no entries
     if (flag == 0){
       myOut = yield_int
       flag = 1
     } else {
       myOut2 = yield_int
       myOut = rbind(myOut, myOut2)
     }
   }
   # clean
  rm(myOut2)
  rm(asappliedOut2)
  rm(asplantedOut2)
  rm(yield_int)
  rm(asapplied_int)
  rm(asplanted_int)
  rm(samps)
 }


 # reassign
 #yield <- myOut
 #nitrogen <- asappliedOut
 #planting <- asplantedOut

 my_list <- list("yield" = myOut, "asapplied" = asappliedOut, "asplanted" = asplantedOut)
 return(my_list)

}

# JPN: download data script
download_workshop_data <- function(reDownloadData = TRUE, reDownloadTrialData = TRUE, downloadScripts = TRUE){
  #reDownloadData = TRUE # only if we wanna re-download the data for any reason
  #reDownloadTrialData = TRUE # this will download "new" trial data that was simulated
  #downloadScripts = TRUE # download lesson catchups

  # URLs of data: original data
  dataURLS = c("https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/asapplied.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/asplanted.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/abline.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/boundary.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/trial.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/yield.gpkg",
	       "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/data/fertilizer_use.csv")

  # URLs of simulated data
  simsURLS = c("https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/asapplied_new.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/asplanted_new.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/trial_new.gpkg",
	       "https://github.com/data-carpentry-for-agriculture/trial-lesson/raw/gh-pages/_episodes_rmd/data/yield_new.gpkg")

  # URLS of catchup scripts
  scriptURLS = c("https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson1.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson2.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson3.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson4.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson5.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson6.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson7.R",
		 "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/lesson-catchup-scripts/lesson8.R")

  # functions.R download
  functionsURL = "https://raw.githubusercontent.com/data-carpentry-for-agriculture/trial-lesson/gh-pages/_episodes_rmd/functions.R"

  # download the data
  if (dir.exists(paste0(getwd(),"/data"))){
    #print('Data directory exists!')
    if (reDownloadData){
      for (i in 1:length(dataURLS)){
      	# get names of files
      	myList = strsplit(dataURLS[i],'/')
      	fname = tail(myList[[1]],1)
      	download.file(dataURLS[i], paste0(getwd(),"/data/",fname), method = "auto", quiet=FALSE, mode="wb")
      }
    }
    # also update new ones
    if (reDownloadTrialData){
      for (i in 1:length(simsURLS)){
      	# get names of files
      	myList = strsplit(simsURLS[i],'/')
      	fname = tail(myList[[1]],1)
      	download.file(simsURLS[i], paste0(getwd(),"/data/",fname), method = "auto", quiet=FALSE, mode="wb")
      }
    }
  } else { # create directory if doesn't exist
    dir.create(paste0(getwd(),"/data"))
    # grab all data
    for (i in 1:length(dataURLS)){
      # get names of files
      myList = strsplit(dataURLS[i],'/')
      fname = tail(myList[[1]],1)
      download.file(dataURLS[i], paste0(getwd(),"/data/",fname), method = "auto", quiet=FALSE, mode="wb")
    }
    if (reDownloadTrialData){
      for (i in 1:length(simsURLS)){
      	# get names of files
      	myList = strsplit(simsURLS[i],'/')
      	fname = tail(myList[[1]],1)
      	download.file(simsURLS[i], paste0(getwd(),"/data/",fname), method = "auto", quiet=FALSE, mode="wb")
      }
    }
  }

  # Also download all catchup scripts
  if (downloadScripts){
    if (!(dir.exists(paste0(getwd(),"/lesson-catchup-scripts")))){ # check if directory exists, make it otherwise
      dir.create(paste0(getwd(),"/lesson-catchup-scripts"))
    }

    for (i in 1:length(scriptURLS)){
      # get names of files
      myList = strsplit(scriptURLS[i],'/')
      fname = tail(myList[[1]],1)
      download.file(scriptURLS[i], paste0(getwd(),"/lesson-catchup-scripts/",fname), method = "auto", quiet=FALSE)
    }
  }

  # finally, download functions.R
  myList = strsplit(functionsURL,'/')
  fname = tail(myList[[1]],1)
  download.file(functionsURL, paste0(getwd(),'/', fname), method = "auto", quiet=FALSE)
}

## JPN: run test script
run_workshop_test <- function(workingDir="default"){
  # copied from before:
  library("rgdal")
  library("plyr")
  library("dplyr")
  library("sp")
  library("sf")
  library("gstat")
  library("tmap")
  library("measurements")
  library("daymetr")
  library("FedData")
  library("lubridate")
  library("raster")
  library("data.table")
  library("broom")
  library("ggplot2")

  #set margins
  par(mar=c(1,1,1,1))

  # create working directory
  if (workingDir != "default"){ # if not default, store there
    print(paste0('You are going to install things in ', workingDir, "/WorkingDir"))
    if (!(dir.exists(paste0(workingDir,"/WorkingDir")))){
      dir.create(paste0(workingDir,"/WorkingDir"))
      # set our working directory to what we have just created
      setwd(paste0(workingDir,"/WorkingDir"))
    }
  } else { # default
    print("Creating file: WorkingDir")
    if (!(dir.exists(paste0(getwd(),"/WorkingDir")))){
      dir.create(paste0(getwd(),"/WorkingDir"))
      # set our working directory to what we have just created
      setwd(paste0(getwd(),"/WorkingDir"))
    }
  }

  # make test images
  if (!(dir.exists(paste0(getwd(),"/test_images")))){ # where to save test images, create dir if not there
    dir.create(paste0(getwd(),"/test_images"))
  }

  # plot #1
  jpeg(paste0(getwd(),"/test_images/plot1_ggplot2.jpg"))
  # test ggplot2
  df = data.frame("x"=c(1,2,3), "y"=c(1,2,3))
  myPlot1 = ggplot(df, aes(x=x, y=y)) + geom_point()
  #plot(myPlot1)
  #pushViewport(myPlot1)
  # grab data, plot with default
  plot(myPlot1$data)
  dev.off()


  # plot #2
  jpeg(paste0(getwd(),"/test_images/plot2_sf.jpg"))
  # test sf
  plot.new()
  nc <- st_read(system.file("shape/nc.shp", package="sf"))
  #pushViewport(plot(nc$geometry) )
  #x = recordPlot(nc$geometry)
  plot(nc$geometry)
  dev.off()

  # text #1
  # test dplyr
  names(nc)[1]
  nc <- dplyr::rename(nc, area = AREA)
  text1 = names(nc)[1]

  # plot #3
  jpeg(paste0(getwd(),"/test_images/plot3_gstat.jpg"))
  # test gstat
  data(meuse)
  coordinates(meuse) = ~x+y
  data(meuse.grid)
  coordinates(meuse.grid) = ~x+y
  gridded(meuse.grid) = TRUE
  #   plot
  lzn.vgm = variogram(log(zinc)~1,data = meuse)
  lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
  myPlot = plot(lzn.vgm, lzn.fit)
  plot(myPlot)
  dev.off()

  # plot #4
  jpeg(paste0(getwd(),"/test_images/plot4_tmap.jpg"))
  # test tmap
  #nc <- st_read(system.file("shape/nc.shp", package="sf"))
  myMap = tm_shape(nc) + tm_polygons('area')
  #plot(myMap)
  tmap_save(myMap, paste0(getwd(),"/test_images/plot4_tmap.jpg"))
  dev.off()

  # text #2
  # test measurements
  text2 = conv_unit(2.54, "cm", "inch")

  # plot #5
  jpeg(paste0(getwd(),"/test_images/plot5_daymetr.jpg"))
  # test daymetr
  df <- download_daymet(site = "Oak Ridge National Laboratories",
			lat = 36.0133,
			lon = -84.2625,
			start = 2000,
			end = 2010,
			internal = TRUE,
			simplify = TRUE)
  plot(df$year, df$value)
  dev.off()

  # plot 6
  jpeg(paste0(getwd(),"/test_images/plot6_ssurgo.jpg"))
  # test FedData
  ssurgo <- get_ssurgo(template=c('NC019'), label='county')
  plotssurgo <- ssurgo$spatial
  plotssurgo <- plotssurgo[4000,]
  plot(plotssurgo)
  dev.off()

  # text #3
  # test lubridate
  text3 = ymd("20110604")

  # plot 7
  # test raster
  jpeg(paste0(getwd(),"/test_images/plot7_raster.jpg"))
  r <- raster(nrows=10, ncols=10)
  r <- setValues(r, 1:ncell(r))
  plot(r)
  dev.off()

  # plot 8
  jpeg(paste0(getwd(),"/test_images/plot8_datatable.jpg"))
  # test data.table
  DT = data.table(
    ID = c("b","b","b","a","a","c"),
    a = 1:6,
    b = 7:12,
    c = 13:18
  )
  #DT
  plot(DT$a, DT$b)
  dev.off()

  # plot 9
  jpeg(paste0(getwd(),"/test_images/plot9_broom.jpg"))
  # test broom
  lmfit <- lm(mpg ~ wt, mtcars)
  myFit = tidy(lmfit)
  plot(c(myFit$estimate, myFit$statistic, myFit$std.error, myFit$p.value))
  dev.off()

  print(' ')
  print(' ')
  print('Three text elements:')
  print(text1)
  print(text2)
  print(text3)
  dev.off()
}

profit_graphs <- function(data, s_ls, n_ls, s_sq, n_sq, Pc, Ps, Pn, other_costs){
  model <- mgcv::gam(yield~
                       s(s, k = 3) +
                       s(n, k = 3),
                     data = data)

  beta <- model$coef
  V_beta <- model$Vp

  data_new <- expand.grid(s_ls, n_ls) %>%
    data.table() %>%
    setnames(names(.),c('s','n'))

  yhat <- predict(model, newdata = data_new, se.fit=TRUE)

  data_pi <- data_new %>%
    .[,y_hat:=yhat$fit] %>%
    .[,y_hat_se:=yhat$se.fit] %>%
    .[,pi_hat:=Pc*y_hat-Ps*s-Pn*n - other_costs] %>%
    .[,pi_hat_se:=Pc*y_hat_se]

  data_seed <- data_pi[n == n_sq,]
  data_nitrogen <- data_pi[s == s_sq,]

  opt_s <- data_seed[,.SD[which.max(pi_hat),]][,s]
  opt_n <- data_nitrogen[,.SD[which.max(pi_hat),]][,n]

  #--- get the model matrix for seed---#
  base_X <- predict(model, newdata = data_seed[s == opt_s,],type='lpmatrix')
  comp_data <- data_seed[s != opt_s,]
  comp_X <- predict(model, newdata = comp_data,type='lpmatrix')

  for (i in 1:(length(s_ls)-1)){
    temp_X <- base_X - comp_X[i,]
    comp_data[i,a:=Pc*(temp_X %*% beta)]
    comp_data[i,b:=Ps*(opt_s-s)]
    comp_data[i,pi_dif:=a-b]
    comp_data[i,pi_dif_se:=temp_X %*% V_beta %*% t(temp_X) %>% sqrt()*Pc]
  }

  profitdiffs <- ggplot(data=comp_data) +
    geom_point(aes(y=-pi_dif,x=factor(s)), color='red') +
    geom_errorbar(aes(ymin=-(pi_dif-1.96*pi_dif_se), ymax=-(pi_dif+1.96*pi_dif_se),x=factor(s)), width=.2,position=position_dodge(.9)) +
    ylab('Profit Diff ($/acre)') +
    xlab('Seed Rate')
  profitdiffs

  # making nitrogen graph
  base_X <- predict(model, newdata=data_nitrogen[n==opt_n,], type='lpmatrix')
  comp_data <- data_nitrogen[n!=opt_n,]
  comp_X <- predict(model, newdata=comp_data, type='lpmatrix')

  for (i in 1:(length(n_ls)-1)){
    temp_X <- base_X - comp_X[i,]
    comp_data[i,a:=Pc*(temp_X %*% beta)]
    comp_data[i,b:=Pn*(opt_n-n)]
    comp_data[i,pi_dif:=a-b]
    comp_data[i,pi_dif_se:=temp_X %*% V_beta %*% t(temp_X) %>% sqrt()*Pc]
  }

  profitdiffn <- ggplot(data=comp_data) +
    geom_point(aes(y=-pi_dif,x=factor(n)), color='red') +
    geom_errorbar(aes(ymin=-(pi_dif-1.96*pi_dif_se), ymax=-(pi_dif+1.96*pi_dif_se),x=factor(n)), width=.2,position=position_dodge(.9)) +
    ylab('Profit Diff ($/acre)') +
    xlab('Nitrogen Rate')
  profitdiffn

  ## now for the profit response curves
  s.seq = seq(min(s_ls), max(s_ls), by = 500)
  n.seq = seq(min(n_ls), max(n_ls), by = 10)

  data_new <- expand.grid(s.seq, n.seq) %>%
    data.table() %>%
    setnames(names(.),c('s','n'))

  yhat <- predict(model, newdata = data_new, se.fit=TRUE)

  data_pi <- data_new %>%
    .[,y_hat:=yhat$fit] %>%
    .[,y_hat_se:=yhat$se.fit] %>%
    .[,pi_hat:=Pc*y_hat-Ps*s-Pn*n - other_costs] %>%
    .[,pi_hat_se:=Pc*y_hat_se]

  data_seed <- data_pi[n == n.seq[(round(length(n.seq)/2))],]
  data_nitrogen <- data_pi[s == s.seq[(round(length(s.seq)/2))],]

  profits <- ggplot() +
    geom_smooth(data=data_seed, method = "loess", aes(y=pi_hat, x=s), size = 0.5,se=FALSE) +
    ylab('Profit ($/acre)') +
    xlab('Seed') +
    theme_grey(base_size = 12)
  profits

  profitn<- ggplot() +
    geom_smooth(data=data_nitrogen, method = "loess", aes(y=pi_hat,x=n), size = 0.5,se=FALSE) +
    ylab('Profit ($/acre)') +
    xlab('Nitrogen') +
    theme_grey(base_size = 12)

  graphs <- list(profitdiffs, profits, profitdiffn, profitn)

  return(graphs)
}

month_prec_graph <- function(prec_merged){
  prec_merged$prec_diff <- prec_merged$prec_month - prec_merged$prec_avg

  monthly_prec <- ggplot(prec_merged) +
    geom_bar(aes(x = month, y = prec_month, fill=prec_diff), stat = 'identity') +
    ggtitle("2018 Monthly Precipitation Compared to Average") +
    labs(y= "Precipitation (in)", x = "Month") +
    guides(fill=guide_legend(title="2018 Minus Average"))
  return(monthly_prec)
}
