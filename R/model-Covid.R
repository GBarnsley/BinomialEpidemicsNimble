#' The epidemic class specific to an SIR model of COVID-19 outbreak in the UK
#' @export
COVIDUKclass <- setClass(
  "COVIDUK",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
#' Function to create an epidemic model of the COVID-19 outbreak in the UK.
#' @param newD Matrix of the national timeseries data, a row per region and time in the columns.
#' @param Pop A vector of region populations, must be same order as newD rows.
#' @param Connectivity A matrix of the transformed distances between region centres.
#' @param ChangePoint Positions in the timeseries that the two lockdowns came into effect.
#' @param TestCapacity Timeseries of the test capacities of entire country.
#' @param StartRegion Index of the starting region.
#' @param TotalInfections The full size of the epidemic, including unobserved infections.
#' @param t.step The time step of the model.
#' @param Frequency TRUE/FALSE value specifying if the model has frequency based transmission.
#' @return An object of COVID class with the compiled model code
#' @export
COVIDModel <- function(newD,
                 Pop,
                 TestCapacity,
                 Connectivity,
                 ChangePoint,
                 StartRegion,
                 t.step = 1,
                 Frequency = TRUE){
  tempCode <- nimbleCode({
    # Model Parameters
    Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
    Alpha ~ dgamma(shape = AlphaShape, rate = AlphaRate)
    Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
    Lockdown[1] ~ dgamma(shape = LockdownShape[1], rate = LockdownRate[1])
    Lockdown[2] ~ dgamma(shape = LockdownShape[2], rate = LockdownRate[2])
    # Likelihood
    for(region in 1:Regions){
      S[region,1] <- Pop[region] - (region == StartRegion)
      I[region,1] <- 0 + (region == StartRegion)
      for (t in 1:TimePeriod){
        newI[region,t] ~ dbinom(size = S[region,t],
                                prob =  probGen(
                                  sum(connectivity[region,1:Regions]*I[1:Regions,t]/((Pop)^Frequency))*
                                    Beta*
                                    Lockdown[1]^((t>=ChangePoint[1] & t<ChangePoint[2])|(t>=ChangePoint[3] & t<ChangePoint[4]))*
                                    Lockdown[2]^((t>=ChangePoint[2] & t<ChangePoint[3])|t>=ChangePoint[4])*
                                    t.step))
        newD[region,t] ~ dbinom(size = I[region,t], prob = probGen(Alpha*t.step*TestCapacity[t]))
        newR[region,t] ~ dbinom(size = I[region,t] - newD[region,t], prob = probGen(Gamma*t.step))
        S[region,t+1] <- S[region,t] - newI[region,t]
        I[region,t+1] <- I[region,t] + newI[region,t] - newD[region,t] - newR[region,t]
      }
    }
  })
  return(COVIDUKclass(
    Model = compileNimble(
      nimbleModel(
        code = tempCode,
        constants = list(TimePeriod = ncol(newD),
                         ChangePoint = ChangePoint,
                         Regions = nrow(newD),
                         t.step = t.step,
                         Connectivity = Connectivity
                         ),
        data = list(newD = newD,
                    Pop = Pop,
                    StartRegion = StartRegion,
                    Frequency = Frequency,
                    TestCapacity = TestCapacity,
                    GammaRate = 1,
                    GammaShape = 1,
                    AlphaRate = 1,
                    AlphaShape = 1,
                    BetaRate = 1,
                    BetaShape = 1,
                    LockdownRate = c(1,1),
                    LockdownShape = c(1,1)
        ),
        inits = list(Beta = 1,
                     Gamma = 1,
                     Alpha = 1,
                     Lockdown = c(0.8,0.9),
                     newI = newD,
                     newR = newD
        ),
        calculate = FALSE
      )
    ),
    MCMC = NA,
    Samples = NA
  )
  )
}
#' Initialization method for the COVID-19 UK SIR model.
#' Sets up provided initial values and calculates the initial
#' newI and newR for those values.
#' @param epiModel An object of the COVIDUK class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return COVIDUK class with the initial values
#' @export
initialValues.COVIDUK <- function(epiModel, hyperParameters){
  #StartingValues
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Alpha <- hyperParameters$`Initial Values`$Alpha
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$Lockdown <- hyperParameters$`Initial Values`$Lockdown
  #Priors
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Shape
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Rate
  epiModel@Model$AlphaShape <- hyperParameters$Priors$Alpha$Shape
  epiModel@Model$AlphaRate <- hyperParameters$Priors$Alpha$Rate
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  epiModel@Model$LockdownRate <- hyperParameters$Priors$Lockdown$Rate
  epiModel@Model$LockdownShape <- hyperParameters$Priors$Lockdown$Shape
  #generating initial newI and newR
  newD <- epiModel@Model$newD
  newR <- round(newD*hyperParameters$Priors$ProportionUndetected)
  #storing first row of dections/removals
  firstRow <- newD[,1] + newR[,1]
  #setting newI as timeseries of detections/removals without first row
  newI <- cbind(newD[,-1],0) + cbind(newR[,-1],0)
  #readding first row
  newI[,1] <- newI[,1] + firstRow
  #removing an infection from start region since we assume it begins with an already exiting infection
  newI[epiModel@Model$StartRegion,newI[epiModel@Model$StartRegion,]!=0][1] <- newI[epiModel@Model$StartRegion,newI[epiModel@Model$StartRegion,]!=0][1] - 1
  epiModel@Model$newR <- newR
  epiModel@Model$newI <- newI
  return(
    epiModel
  )
}
#' Method to build an MCMC for the COVIDUK class.
#' Applies a block RW to the beta, alpha, gamma and lockdown parameters and a NpmDelta algorithm
#' on every region for newI and newR.
#' @param epiModel An object of the COVIDUK class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return a complied MCMC
#' @export
buildMCMCInternal.COVIDUK <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  for(i in 1:nrow(epiModel@Model$newD)){
    output$addSampler(target = paste0("newI[",i,",]"),
                      type = sampler,
                      control = hyperParameters$`N+-Delta`
                      )
    output$addSampler(target = paste0("newR[",i,",]"),
                      type = sampler,
                      control = hyperParameters$`N+-Delta`
                      )
  }
  output$addSampler(target = "Alpha",
                    type = "RW")
  output$addSampler(target = "Beta",
                    type = "RW")
  output$addSampler(target = "Lockdown[1]",
                    type = "RW")
  output$addSampler(target = "Lockdown[2]",
                    type = "RW")
  output$addSampler(target = "Gamma",
                    type = "RW")
  output$addMonitors(c('Alpha','Beta', 'Gamma','Lockdown[1]','Lockdown[2]','newI','newR'))
  output <- buildMCMC(output)
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}
#' A ultility function to import national timeseries data for use in COVIDUK model.
#' @param population a data frame of populations, must include region codes under "code"
#' and populations under "pop".
#' @param positions a data frame of the centres of the regions, must include region codes,
#' under "code", y-coord as "y" and x-coord as "x".
#' @param StartRegion Region Code where the first case occured
#' @param Lockdown1Date Date first lockdown occured
#' @param Lockdown2Date Date second lockdown occured
#' @param startDate the date of the first infection
#' @return a list of timeseries
#' @export
ImportCOVIDUKTimeSeries <- function(population,
                                    positions,
                                    StartRegion = "E12000001",
                                    Lockdown1Date = lubridate::dmy("23/03/2020"),
                                    Lockdown2Date = lubridate::dmy("13/05/2020"),
                                    Lockdown1Date2 = lubridate::dmy("01/11/2020"),#CHECK THIS!
                                    Lockdown2Date2 = lubridate::dmy("01/12/2020"),#CHECK THIS!
                                    VaccinationStartDate = lubridate::dmy("08/12/2020"),
                                    startDate = lubridate::dmy("20/01/20")){
  ###Loading Case Data
  ##Loading England Region Cases:
  EnglandRawCases <- loadCOVIDAPIdate("region","newCasesBySpecimenDate","newCases")
  ##Loading Scotland and Wales
  NationRawCases <- loadCOVIDAPIdate("nation","newCasesBySpecimenDate","newCases")
  ##Combing Cases
  RawCases <- rbind(EnglandRawCases, NationRawCases[NationRawCases$code %in% c("S92000003","W92000004"),])
  ##storing region codes
  regionCodes <- sort(unique(RawCases$code))
  ##storing Regions (the constant)
  Regions <- length(regionCodes)
  ##Storing max date in timeseries (used later in determining overall length), we subtract 7 because of delays in reporting
  CasesMaxDate <- max(RawCases$date)-7
  ##Sorting into matrix (newD)
  newD <- matrix(0,nrow = length(regionCodes), ncol = max(RawCases$date)-startDate)
  for(row in 1:nrow(RawCases)){
    newD[which(RawCases$code[row] == regionCodes), RawCases$date[row] - startDate] <- RawCases$newCases[row]
  }
  ###Setting Up Connectivity
  Connectivity <- matrix(NA, nrow = Regions, ncol = Regions)
  for(regionIndex1 in 1:Regions){
    for(regionIndex2 in 1:Regions){
      xy1 <- c(positions$x[positions$code == regionCodes[regionIndex1]], 
               positions$y[positions$code == regionCodes[regionIndex1]])
      xy2 <- c(positions$x[positions$code == regionCodes[regionIndex2]], 
               positions$y[positions$code == regionCodes[regionIndex2]])
      Connectivity[regionIndex1,regionIndex2] <- sqrt(sum((xy1 - xy2)^2))
    }
  }
  #standardising
  Connectivity <- Connectivity*(1/max(Connectivity))
  #inverting
  Connectivity <- 1-Connectivity
  ###Setting up population counts
  ##just put populations in correct order
  pop <- rep(NA, Regions)
  for(regionIndex in 1:Regions){
    pop[regionIndex] <- population$pop[population$code == regionCodes[regionIndex]]
  }
  ###Loading Test Capacity Data
  RawTestCapacity <- loadCOVIDAPIdate("overview","plannedCapacityByPublishDate","test")
  ##Storing Maximum date
  TestMaxDate <- max(RawTestCapacity$date)
  #flipping to correct order and imputing capacity from start date to first set of date, just assume linear increase from 0
  TestCapacity <- c(seq(0,RawTestCapacity$test[which.min(RawTestCapacity$date)],length.out=min(RawTestCapacity$date)-startDate),
                    rev(RawTestCapacity$test))
  ###Setting up ChangePoints
  ChangePoint <- c(Lockdown1Date, Lockdown2Date, Lockdown1Date2, Lockdown2Date2) - startDate
  ###Setting up StartRegion
  StartRegion <- which(regionCodes == StartRegion)
  ###Determining TimePeriod
  TimePeriod <- min(c(CasesMaxDate,TestMaxDate,VaccinationStartDate)) -startDate
  ##Trimming timeseries
  newD <- newD[,1:TimePeriod]
  TestCapacity <- TestCapacity[1:TimePeriod]
  ###Setting up the output lists
  Output <- list(
    newD = newD,
    Pop = pop,
    TestCapacity = TestCapacity,
    Connectivity = Connectivity,
    ChangePoint = ChangePoint,
    StartRegion = StartRegion
  )
  return(Output)
}
#' Internal function that loads data from the government api
#' @param regionType a string specifiying the type of region wanted e.g. "nation" "region"
#' @param output a string specifying the output wanted e.g. "newCasesBySpecimenDate" "plannedCapacityByPublishDate"
#' @param outputName the name of the output in the return data frame.
#' @return a data frame of the dates, codes and output
#' @export
loadCOVIDAPIdate <- function(regionType, output, outputName){
  endpoint <- paste0('https://api.coronavirus.data.gov.uk/v1/data?filters=areaType=',
    regionType,'&structure={"code":"areaCode","date":"date","',
    outputName,'":"',
    output,'"}')
  httr::GET(url = endpoint, httr::timeout(10)) -> response
  if (response$status_code >= 400) {
    err_msg = httr::http_status(response)
    stop(err_msg)
  }
  json_text <- httr::content(response, "text")
  output <- jsonlite::fromJSON(json_text)$data
  output$date <- lubridate::ymd(output$date)
  return(output)
}
