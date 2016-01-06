## DPLYR ##
## dplyr ##

install.packages("dplyr")
library(dplyr)

install.packages("nycflights13")
library(nycflights13)

dim(flights)

head(flights)


## ALTERNATIVELY ##

install.packages("dplyr")
library(dplyr)

setwd('/Users/clinton/Documents/DataScience/Day6-R')

flights <- read.csv("input_files/nycflights13_flights.csv")

dim(flights)

head(flights)



## CONVERT TO DPLYR DATA FRAME ##

flights <- tbl_df(flights)

glimpse(flights)

flights



## SELECTING COLUMNS ##

select(flights, carrier)

select(flights, carrier, flight)

select(flights, year:day)

select(flights, -distance)

select(flights, -(year:dest))

select(flights, contains("time"))

select(flights, starts_with("air"))

select(flights, ends_with("delay"))

select(flights, matches("^arr_.*"))

select(flights, one_of("year", "month", "day", "carrier"))

select(flights, 2:4)


## RENAMING COLUMNS ##

rename(flights, tail_num=tailnum)

rename(flights, destination=dest)


## ADDING NEW COLUMNS ##

mutate(flights,
  gain = arr_delay - dep_delay,
  speed = distance / air_time * 60)

mutate(flights,
  gain = arr_delay - dep_delay,
  gain_per_hour = gain / (air_time / 60)
)



## FILTERING ROWS ##

## BY POSITION ##

slice(flights, 1:10)

slice(flights, 1:10)[, 2:4]


## BY CONDITION(S) ##

filter(flights, month==1, day==1)

filter(flights, month==1 | month==2)

filter(flights, (month==1 | month==2) & carrier == "AA")

filter(flights, carrier == "UA") %>% select(origin, dest)

filter(flights, dep_delay > 60) %>% select(year:day, carrier, dep_delay)


## ARRANGING ROWS ##

arrange(flights, year, month, day) %>% select(year, month, day)

arrange(flights, desc(arr_delay))[ , c(8,9,10)]

arrange(flights, desc(distance)) %>% select(origin, dest, distance)


## DISTINCT ROWS ##

distinct(select(flights, tailnum))

distinct(select(flights, origin, dest))


## RANDOM SAMPLE OF ROWS ##

sample_n(flights, 10)

sample_n(flights, 10, replace=TRUE)

sample_frac(flights, 0.01)

sample_frac(flights, 0.01, replace=TRUE)



## SUMMARIZE WITH SUMMARISE ##

summarise(flights,
  delay = mean(dep_delay, na.rm = TRUE))



## GROUP BY ##

destinations <- group_by(flights, dest)
summarise(destinations,
  planes = n_distinct(tailnum),
  flights = n()
)


by_tailnum <- group_by(flights, tailnum)
delay <- summarise(by_tailnum,
  count = n(),
  dist = mean(distance, na.rm = TRUE),
  delay = mean(arr_delay, na.rm = TRUE))
delay <- filter(delay, count > 20, dist < 2000)

install.packages("ggplot2")
library(ggplot2)
ggplot(delay, aes(dist, delay)) +
  geom_point(aes(size = count), alpha = 1/2) +
  geom_smooth() +
  scale_size_area()


daily <- group_by(flights, year, month, day)
(per_day   <- summarise(daily, flights = n()))

(per_month <- summarise(per_day, flights = sum(flights)))

(per_year  <- summarise(per_month, flights = sum(flights)))



## CHAINING ##

Use chaining to go from this step-by-step code:

a1 <- group_by(flights, year, month, day)
a2 <- select(a1, arr_delay, dep_delay)
a3 <- summarise(a2,
  arr = mean(arr_delay, na.rm = TRUE),
  dep = mean(dep_delay, na.rm = TRUE))
a4 <- filter(a3, arr > 30 | dep > 30)

To chained code that doesnt require you to save the intermediate results:

flights %>%
  group_by(year, month, day) %>%
  select(arr_delay, dep_delay) %>%
  summarise(
    arr = mean(arr_delay, na.rm = TRUE),
    dep = mean(dep_delay, na.rm = TRUE)
  ) %>%
  filter(arr > 30 | dep > 30)






## DATA.TABLE ##
## data.table ##

install.packages("data.table")
library(data.table)

install.packages("nycflights13")
library(nycflights13)

dim(flights)

head(flights)


## ALTERNATIVELY ##

install.packages("data.table")
library(data.table)

setwd("~/Documents/Data-Science-with-Unix-SQL-R-Python/Day6-R")

flights <- fread("input_files/nycflights13_flights.csv")

dim(flights)

head(flights)



## SELECTING COLUMNS ##

flights[,carrier]

flights[,.(carrier,flight)]   # OR    flights[,list(carrier,flight)]

#select(flights, year:day)
flights[,.(year, month, day)]

#select(flights, -distance)

#select(flights, -(year:dest))

#select(flights, contains("time"))

#select(flights, starts_with("air"))

#select(flights, ends_with("delay"))

#select(flights, matches("^arr_.*"))

#select(flights, one_of("year", "month", "day", "carrier"))

#select(flights, 2:4)
flights[,.(year, month, day)]


## RENAMING COLUMNS ##

setnames(flights, "tailnum", "tail_num")
setnames(flights, "tailnum", "tail_num")[]

setnames(flights, c("tail_num","dest"), c("tail_number","destination"))
setnames(flights, c("tail_num","dest"), c("tail_number","destination"))[]


## ADDING NEW COLUMNS ##

flights[, gain := arr_delay - dep_delay]

flights[, c("gain","speed") := list(arr_delay - dep_delay, 
									distance / air_time * 60)]

flights[, c("gain","gain_per_hour") := list(arr_delay - dep_delay, 
											gain / (air_time / 60))]

flights[, V1 := NULL]    # DELETE A COLUMN



## FILTERING ROWS ##

## BY POSITION ##

flights[1:10]        # OR    flights[1:10, ]

flights[1:10, .(year, month, day)]


## BY CONDITION(S) ##

flights[month==1 & day==1]
# OR    setkey(flights, month, day)    flights[.(1,1)]

flights[month==1 | month==2]
# OR    setkey(flights, month)    flights[.(1|2)]

flights[(month==1 | month==2) & carrier=="AA"]
# OR    setkey(flights, month, carrier)    flights[.(1|2, "AA")]

flights[carrier=="UA"][,.(origin, destination)]
# OR    setkey(flights, carrier)    flights["UA"][,.(origin, destination)]

flights[dep_delay > 60][,.(year, month, day, carrier, dep_delay)]


## ARRANGING ROWS ##

flights[order(year,month,day)][,.(year,month,day)]

flights[order(-arr_delay)][,.(carrier, flight, origin, destination, arr_delay)]

flights[order(-distance)][,.(origin, destination, distance)]


## DISTINCT ROWS ##

unique(flights[,.(tail_number)])        # OR   unique(flights[,tail_number])
unique(flights[,.(origin,destination)])    # unique(flights[,c(origin,destination)])
unique(flights[,.(origin,destination)])[order(origin)]
unique(flights[,.(origin,destination)])[order(origin,destination)]


## RANDOM SAMPLE OF ROWS ##

flights[sample(.N, 10)]

flights[sample(.N, 10, replace=TRUE)]
flights[,.(origin, destination),by=carrier][sample(.N, 10, replace=FALSE)]

sample_frac(flights, size=0.01)

sample_frac(flights, size=0.01, replace=TRUE)


## SUMMARIZE ##

flights[,.(delay=mean(dep_delay,na.rm=TRUE))]



## GROUP BY ##

flights[,.(planes = length(unique(tail_number)), flights = .N), by=destination][order(destination)]


delay <- flights[,.(count=.N, dist = mean(distance, na.rm=TRUE),
  delay = mean(arr_delay, na.rm=TRUE)), by=tail_number][count > 20 & dist < 2000][order(tail_number)]

install.packages("ggplot2")
library(ggplot2)
ggplot(delay, aes(dist, delay)) +
  geom_point(aes(size = count), alpha = 1/2) +
  geom_smooth() +
  scale_size_area()


flights[,.(flights=.N),by=.(year,month,day)][order(year,month,day)]

flights[,.(flights=.N),by=.(year,month)][order(year,month)]

flights[,.(flights=.N),by=year][order(year)]



## CHAINING ##

Use chaining to go from this step-by-step code:

a1 <- group_by(flights, year, month, day)
a2 <- select(a1, arr_delay, dep_delay)
a3 <- summarise(a2,
  arr = mean(arr_delay, na.rm = TRUE),
  dep = mean(dep_delay, na.rm = TRUE))
a4 <- filter(a3, arr > 30 | dep > 30)

To chained code that doesnt require you to save the intermediate results:

flights[,.(arr = mean(arr_delay, na.rm = TRUE),
dep = mean(dep_delay, na.rm = TRUE)),by=.(year,month,day)][arr > 30 | dep > 30][order(year,month,day)]






## GRAPHING AND PLOTTING ##

## GGPLOT2 ##


## SAVING A GRAPH ##

## FROM SCRIPT OR FUNCTION ##

pdf("plots.pdf")
print(qplot(...))
dev.off()

## FROM SCREEN ##

ggsave("plots.png", width=5, height=5, dpi=100)



## HISTOGRAMS ##

summary(flights$distance)
## Basic histogram from the vector "distance". Each bin is 200 wide.
ggplot(flights, aes(x=distance)) + geom_histogram(binwidth=200)

# Draw with black outline, white fill
ggplot(flights, aes(x=distance)) +
    geom_histogram(binwidth=200, colour="black", fill="white")

# Add a line for the mean
ggplot(flights, aes(x=distance)) +
    geom_histogram(binwidth=200, colour="black", fill="white") +
    geom_vline(aes(xintercept=mean(distance, na.rm=T)),
               color="red", linetype="dashed", size=1)

library(data.table)
ua_aa <- data.table(flights[flights$carrier %in% c('UA', 'AA'), ])
# Overlaid histograms
ggplot(ua_aa, aes(x=distance, fill=carrier)) +
    geom_histogram(binwidth=200, alpha=.5, position="identity")

ua_aa_means <- ua_aa[,.(distance_mean = mean(distance, na.rm = TRUE)),by=carrier][order(distance_mean)]
# Overlaid histograms with means
ggplot(ua_aa, aes(x=distance, fill=carrier)) +
    geom_histogram(binwidth=200, alpha=.5, position="identity") +
    geom_vline(data=ua_aa_means, aes(xintercept=distance_mean,  colour=carrier),
               linetype="dashed", size=1)		

# Using facets, with mean lines
ggplot(ua_aa, aes(x=distance)) + 
			geom_histogram(binwidth=200, colour="black", fill="white") +
			facet_grid(carrier ~ .) + 
			geom_vline(data=ua_aa_means, aes(xintercept=distance_mean), 
						linetype="dashed", size=1, colour="red")


	
## DENSITY PLOTS ##

# Kernel Density Plot
ggplot(flights, aes(x=distance)) + geom_density()



## BAR PLOTS ##

ggplot(data=ua_aa, aes(x=carrier, y=distance, fill=carrier)) +
	geom_bar(stat="identity")

# Add title, narrower bars, fill color, and change axis labels
ggplot(data=ua_aa, aes(x=carrier, y=distance, fill=carrier)) + 
    geom_bar(width=.8, stat="identity") + 
    guides(fill=FALSE) +
    xlab("Carrier") + ylab("Total distance") +
    ggtitle("Total Distance by Carrier")

# Bar graph of counts
ggplot(data=ua_aa, aes(x=carrier)) +
    geom_bar(stat="count")



## LINE CHARTS ##

daily_flights <- flights[carrier %in% c('UA', 'AA'), ][,.(flights=.N),by=.(year,month,day,carrier)][order(year,month,day,carrier)]

daily_flights[, date := as.POSIXct(paste(year,month,day), format="%Y %m %d")]

ggplot(data=daily_flights, aes(x=date, y=flights, group=carrier)) +
    geom_line()

# Change color of both line and points
# Change line type and point type, and use thicker line and larger points
# Change points to circles with white fill
ggplot(data=daily_flights, aes(x=date, y=flights, group=carrier)) + 
    geom_line(colour="red", linetype="solid", size=0.5) + 
    geom_point(colour="red", size=1, shape=21, fill="white")



## BOX PLOTS ##

# A basic box plot
ggplot(ua_aa, aes(x=carrier, y=distance)) + geom_boxplot()

# A basic box with the conditions colored
ggplot(ua_aa, aes(x=carrier, y=distance, fill=carrier)) + geom_boxplot()

# The above adds a redundant legend. With the legend removed:
ggplot(ua_aa, aes(x=carrier, y=distance, fill=carrier)) + geom_boxplot() + guides(fill=FALSE)

# With flipped axes
ggplot(ua_aa, aes(x=carrier, y=distance, fill=carrier)) + geom_boxplot() + guides(fill=FALSE) + coord_flip()



## SCATTER PLOTS ##

flights_sample <- sample_frac(flights, size=0.01, replace=FALSE)

ggplot(flights_sample, aes(x=distance, y=air_time)) +
    geom_point(shape=1, alpha=0.10, na.rm=TRUE)      # Use hollow circles

# Add linear regression line, by default includes 95% confidence region
ggplot(flights_sample, aes(x=distance, y=air_time)) +
    geom_point(shape=1, alpha=0.10, na.rm=TRUE) +
    geom_smooth(method=lm, colour="red", linetype="solid", na.rm=TRUE)

# Don't add shaded confidence region
ggplot(flights_sample, aes(x=distance, y=air_time)) +
    geom_point(shape=1, alpha=0.10, na.rm=TRUE) +
    geom_smooth(method=lm, colour="red", linetype="solid", na.rm=TRUE, 
				se=FALSE)

# Add a loess smoothed fit curve with confidence region
ggplot(flights_sample, aes(x=distance, y=air_time)) +
    geom_point(shape=1, alpha=0.10, na.rm=TRUE) +
    geom_smooth(method=loess, 
				color="blue", linetype="solid", size=1.5, na.rm=TRUE, se=FALSE)
