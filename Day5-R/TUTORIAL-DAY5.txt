## R TUTORIAL ##

## DATA TYPES ##

## VECTORS ##

vector_numeric <- c(1,2,5.3,6,-2,4)
vector_character <- c("one","two","three")
vector_logical <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE)

vector_numeric[c(1,3)] # 1st and 3rd elements of vector


## LISTS ##

# a list with 4 components:
# a string, a numeric vector, a matrix, and a scaler 
list_w <- list(name="Fred", my_numeric_vector=vector_numeric, my_matrix=y, age=15.3)

# a list containing two lists 
list_v <- c(list_w,list_w)

list_v[4]
list_v[[4]] # 4th component of the list_v
list_w[["my_numeric_vector"]] # component named my_numeric_vector in list


## FACTORS ##

# variable gender with 20 "male" entries and 30 "female" entries 
gender <- c(rep("male",20), rep("female", 30)) 
gender <- factor(gender)
# stores gender as 20 1s and 30 2s and associates
# 1=female, 2=male internally (alphabetically)
# R now treats gender as a nominal variable 
summary(gender)
levels(gender)

v1 <- c(1,1,3,2)
e <- c("up", "down", "flat", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
mydata <- data.frame(v1,e,f)
class(mydata$v1)
# variable v1 is coded 1, 2 or 3
# we want to attach value labels 1=red, 2=blue, 3=green
mydata$v1 <- factor(mydata$v1,
						levels = c(1,2,3),
						labels = c("red", "blue", "green"))
class(mydata$v1)
levels(mydata$v1)

# variable rating coded as "large", "medium", "small'
rating <- c("large", "small", "medium", "medium", "medium", "small", "large", "large", "large")
rating <- ordered(rating)
rating

# recodes rating to 1,2,3 and associates
# 1=large, 2=medium, 3=small internally
# R now treats rating as ordinal

# variable y is coded 1, 3 or 5 
# we want to attach value labels 1=Low, 3=Medium, 5=High
y <- c(1,3,3,5)
e <- c("up", "down", "flat", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
mydata <- data.frame(y,e,f)
mydata$y <- ordered(mydata$y,
						levels = c(1, 3, 5),
						labels = c("Low", "Medium", "High"))
class(mydata$y)
mydata$y


## MATRICES ##

# generates 5 x 4 numeric matrix 
y<-matrix(1:20, nrow=5, ncol=4)
y

# another example
cells <- c(1,26,24,68)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2") 
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=TRUE,
  							dimnames=list(rnames, cnames))
mymatrix

mymatrix[,2] # 2nd column of matrix
mymatrix[1,] # 3rd row of matrix 
mymatrix[1:2, 1:2] # rows 1 to 2 and columns 1 to 2


## DATA FRAMES ##

d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
my_data <- data.frame(d,e,f)
names(my_data) <- c("ID","Color","Passed") # variable names
names(my_data)
my_data

my_data[2:3] # columns 2 to 3 of data frame
my_data[c("ID","Passed")] # columns ID and Passed from data frame
my_data$Color # variable Color in the data frame


## USEFUL FUNCTIONS ##

length(object) # number of elements or components
str(object)    # structure of an object 
class(object)  # class or type of an object
names(object)  # names

c(object,object,...)       # combine objects into a vector
cbind(object, object, ...) # combine objects as columns
rbind(object, object, ...) # combine objects as rows 

object     # prints the object

ls()       # list current objects
rm(object) # delete an object

newobject <- edit(object) # edit copy and save as newobject 

seq(from, to, by)  # generate a sequence
indices <- seq(1,10,2)  # indices is c(1, 3, 5, 7, 9)
indices

rep(x, ntimes)	repeat x n times
y <- rep(1:3, 2)  # y is c(1, 2, 3, 1, 2, 3)
y

cut(x, n)	divide continuous variable into factor with n levels 
x <- c(1,1,1,1,2,3,2,23,4,5,6,7,8,7,6,5,6,7,6,7,6,7,6,7,8,9,8,7,6,5,4,9)
y <- cut(x, 5)
y

install.packages("ggplot2")   # install a package
library(ggplot2)     # use the package
require(ggplot2)



## IMPORTING / EXPORTING DATA ##

## CSV FILES ##

# first row contains variable names, comma is separator 
# assign the variable id to row names
# note the / instead of \ on mswindows systems
mydata <- read.table("c:/mydata.csv", header=TRUE, sep=",", row.names="id")

setwd("~/Documents/Data-Science-with-Unix-SQL-R-Python/Day5-R")
air_quality <- read.table("input_files/airquality.csv", header=TRUE, sep=",", row.names="observation")
head(air_quality)


write.table(air_quality, "output_files/airquality.tsv", sep="\t")


## EXCEL FILES ##

# read in the first worksheet from the workbook myexcel.xlsx
# first row contains variable names
install.packages("xlsx")
library(xlsx)
mydata <- read.xlsx("c:/myexcel.xlsx", 1)
# read in the worksheet named mysheet
mydata <- read.xlsx("c:/myexcel.xlsx", sheetName = "mysheet")

setwd("~/Documents/Data-Science-with-Unix-SQL-R-Python/Day5-R")
air_quality <- read.xlsx("input_files/airquality.xlsx", 1)
head(air_quality)


write.xlsx(air_quality, "output_files/airquality.xlsx")


## DATABASES ##

Packages: RODBC, RMySQL, ROracle, RJDBC, RSQLite, sqldf

install.packages("RSQLite")
library(RSQLite)
db <- dbConnect(SQLite(), "input_files/DataScience.db")
summary(db)
dbListTables(db)
dbListFields(db, "iris")
iris <- dbReadTable(db, "iris", row.names=NULL, check.names=TRUE, select.cols="sepal_length,sepal_width,petal_length,petal_width,species")
setosa <- iris[iris$species=='setosa', ]
dbWriteTable(db, "Rsetosa", setosa[1:10, ])
dbReadTable(db, "Rsetosa")


## VIEWING DATA ##

# list objects in the working environment
ls()

# list the variables in iris
iris <- read.table("input_files/iris.csv", header=TRUE, sep=',')
iris <- iris[ ,c("sepal_length", "sepal_width", "petal_length", "petal_width", "species")]

names(iris)

# list the structure of iris
str(iris)

# list levels of factor species in iris
iris$species <- factor(iris$species)
levels(iris$species)

# dimensions of an object
dim(iris)

# class of an object (numeric, matrix, data frame, etc)
class(iris)
class(iris$sepal_width)
class(iris$species)

# print iris data
iris

# print first 10 rows of iris
head(iris, n=10)

# print last 5 rows of iris
tail(iris, n=5)



## MISSING DATA ##

## TESTING FOR ##

is.na(x) # returns TRUE of x is missing
y <- c(1,2,3,NA)
is.na(y) # returns a vector (F F F T)


## RECODING ##

# recode 99 to missing for variable v1
# select rows where v1 is 99 and recode column v1 
mydata$v1[mydata$v1==99] <- NA

iris$sepal_length[iris$sepal_length > 7.0] <- NA
iris[is.na(iris$sepal_length), ]


## EXCLUDING ##

x <- c(1,2,NA,3)
mean(x) # returns NA
mean(x, na.rm=TRUE) # returns 2

mean(iris$sepal_length)
mean(iris$sepal_length, na.rm=TRUE)

# list rows of data that have missing values 
iris[!complete.cases(iris),]

# create new dataset without missing data 
iris_no_missing <- na.omit(iris)
dim(iris_no_missing)


## DATES ##

Dates are represented as the number of days since 1970-01-01, with negative values for earlier dates.

# use as.Date( ) to convert strings to dates 
my.dates <- as.Date(c("2015-06-22", "2015-02-13"))

# number of days between 6/22/15 and 2/13/15 
days <- my.dates[1] - my.dates[2]
days

Sys.Date() # returns today's date. 
date() # returns the current date and time.

# print today's date
today <- Sys.Date()
format(today, format="%B %d %Y")


## DATE CONVERSION ##

## CHARACTER TO DATE ##

# convert date info in format 'mm/dd/yyyy'
strDates <- c("01/05/2015", "08/16/2015")
dates <- as.Date(strDates, "%m/%d/%Y")
dates

my_dates <- as.Date(c("2015-12-22", "2015-10-13"))
my_dates

## DATE TO CHARACTER ##

# convert dates to character data
strDates <- as.character(my_dates)
strDates



## CREATING NEW VARIABLES ##

# Three examples for doing the same computations
diamonds <- read.table("input_files/diamonds.csv", header=TRUE, sep=',')
diamonds <- diamonds[ ,c("carat","cut","color","clarity","depth","table","price","x","y","z")]
head(diamonds)

diamonds$sum <- diamonds$x + diamonds$y + diamonds$z
diamonds$mean <- (diamonds$x + diamonds$y + diamonds$z) / 3

rm(x)
rm(y)
rm(z)
attach(diamonds)
diamonds$sum <- x + y + z
diamonds$mean <- (x + y + z) / 3
detach(diamonds)

diamonds$sum <- NA
diamonds$mean <- NA

diamonds <- transform(diamonds,
sum = x + y + z,
mean = (x + y + z) / 3
)

head(diamonds)


## RECODING VARIABLES ##

id <- seq(1,50)
gender <- c(rep("male",20), rep("female", 30))
age <- sample(18:90, size=50, replace=TRUE)
my.data <- data.frame(id, gender, age)
names(my.data) <- c("ID","Gender","Age")
my.data

# create 2 age categories 
my.data$AgeCategory <- ifelse(my.data$Age > 70, c("Older"), c("Younger")) 

# another example: create 3 age categories 
attach(my.data)
my.data$AgeCategories[Age > 75] <- "Elder"
my.data$AgeCategories[Age > 45 & Age <= 75] <- "Middle Aged"
my.data$AgeCategories[Age <= 45] <- "Young"
detach(my.data)
head(my.data)


## RENAMING VARIABLES ##

# rename programmatically
install.packages("reshape")
library(reshape)
my.data <- rename(my.data, c(AgeCategory="OlderYounger"))
names(my.data)

# you can re-enter all the variable names in order
# changing the ones you need to change.the limitation
# is that you need to enter all of them!
names(my.data) <- c("ID","Gender","Age", "OlderYounger", "ThreeAgeCategories")



## ARITHMETIC OPERATORS ##

Operator	Description
+	addition
-	subtraction
*	multiplication
/	division
^ or **	exponentiation
x %% y	modulus (x mod y) 5%%2 is 1
x %/% y	integer division 5%/%2 is 2


## LOGICAL OPERATORS ##

Operator	Description
<	less than
<=	less than or equal to
>	greater than
>=	greater than or equal to
==	exactly equal to
!=	not equal to
!x	Not x
x | y	x OR y
x & y	x AND y
isTRUE(x)	test if X is TRUE



## NUMERIC FUNCTIONS ##

Function	Description
abs(x)	absolute value
sqrt(x)	square root
ceiling(x)	ceiling(3.475) is 4
floor(x)	floor(3.475) is 3
trunc(x)	trunc(5.99) is 5
round(x, digits=n)	round(3.475, digits=2) is 3.48
signif(x, digits=n)	signif(3.475, digits=2) is 3.5
cos(x), sin(x), tan(x)	also acos(x), cosh(x), acosh(x), etc.
log(x)	natural logarithm
log10(x)	common logarithm
exp(x)	e^x


## CHARACTER FUNCTIONS ##

Function	Description
substr(x, start=n1, stop=n2)	Extract or replace substrings in a character vector.
x <- "abcdef" 
substr(x, 2, 4) # is "bcd" 

grep(pattern, x , ignore.case=FALSE, fixed=FALSE)	Search for pattern in x. If fixed=FALSE then pattern is a regular expression. If fixed=TRUE then pattern is a text string. Returns matching indices.
grep("A", c("b","A","c"), fixed=TRUE) # returns 2

sub(pattern, replacement, x, ignore.case =FALSE, fixed=FALSE)	Find pattern in x and replace with replacement text. If fixed=FALSE then pattern is a regular expression.
If fixed = T then pattern is a text string. 
sub("\\s",".","Hello There") # returns "Hello.There"

strsplit(x, split)	Split the elements of character vector x at split. 
strsplit("abc", "") # returns 3 element vector "a","b","c"

paste(..., sep="")	Concatenate strings after using sep string to seperate them.
paste("x",1:3,sep="") # returns c("x1","x2" "x3")
paste("x",1:3,sep="M") # returns c("xM1","xM2" "xM3")
paste("Today is:", date())

x <- "Hello, world!"
toupper(x)	# Uppercase

tolower(x)	# Lowercase


## STATISTICAL FUNCTIONS ##

Function	Description
dnorm(x)	normal density function (by default m=0 sd=1)
# plot standard normal curve
x <- pretty(c(-3,3), 30)
y <- dnorm(x)
plot(x, y, type='l', xlab='x values', ylab='density', yaxs='i')

pnorm(x)	cumulative normal probability for q 
(area under the normal curve to the left of q)
pnorm(1.96) # is 0.975

qnorm(p)	normal quantile. 
value at the p percentile of normal distribution 
qnorm(.9) # is 1.28  90th percentile

rnorm(n, m=0,sd=1)	n random normal variables with mean m and standard deviation sd. 
#50 random normal variables with mean=50, sd=10
x <- rnorm(50, m=50, sd=10)

dbinom(x, size, prob)
pbinom(q, size, prob)
qbinom(p, size, prob)
rbinom(n, size, prob)	binomial distribution where size is the sample size 
and prob is the probability of a heads (pi) 
# prob of 0 to 5 heads of fair coin out of 10 flips
dbinom(0:5, 10, .5) 
# prob of 5 or less heads of fair coin out of 10 flips
pbinom(5, 10, .5)

dpois(x, lamda)
ppois(q, lamda)
qpois(p, lamda)
rpois(n, lamda)	poisson distribution with m=std=lamda
#probability of 0,1, or 2 events with lamda=4
dpois(0:2, 4)
# probability of at least 3 events with lamda=4 
1- ppois(2,4)

dunif(x, min=0, max=1)
punif(q, min=0, max=1)
qunif(p, min=0, max=1)
runif(n, min=0, max=1)	uniform distribution, follows the same pattern 
as the normal distribution above. 
#10 uniform random variates
x <- runif(10)
x
punif(0.73)


Function	Description
mean(x, trim=0, na.rm=FALSE)	mean of object x
# trimmed mean, removing any missing values and 
# 5 percent of highest and lowest scores 
mx <- mean(x, trim=.05, na.rm=TRUE)

sd(x)	standard deviation of object(x).
var(x) for variance
mad(x) for median absolute deviation

median(x)	median

quantile(x, probs)	quantiles where x is the numeric vector whose quantiles are desired and probs is a numeric vector with probabilities in [0,1].
# 30th and 84th percentiles of x
y <- quantile(x, c(.3,.84))
y

range(x)	# range

sum(x)	# sum

diff(x, lag=1)	# lagged differences, with lag indicating which lag to use

min(x)	# minimum

max(x)	# maximum

x_scaled <- scale(x, center=TRUE, scale=TRUE)	# column center or standardize a matrix.
x_scaled


## CONTROL STRUCTURES ##

## IF-ELSE ##

if (cond) expr
if (cond) expr1 else expr2


## FOR ##

for (var in seq) expr

## FOR, IF-ELSE EXAMPLES ##
positives1 <- numeric(); for (element in x_scaled) if (element > 0) positives1 <- c(positives1, element); positives1

positives2 <- numeric(); for (i in 1:length(x_scaled)) if (x_scaled[i] > 0) positives2 <- c(positives2, x_scaled[i]); positives2

positives3 <- numeric(); for (i in 1:length(x_scaled)) if (x_scaled[i] > 0) positives3[i] <- x_scaled[i]; positives3[!is.na(positives3)]


## WHILE ##

while (cond) expr

n <- 1; while (n < 6) { print(n); n = n+1 }


## SWITCH ##

switch(expr, ...)

set.seed(11)
y <- rnorm(5)
x <- "sd"
z <- switch(x,"mean"=mean(y),"median"=median(y),"variance"=var(y),"sd"=sd(y))
z
x <- "median"
z <- switch(x,"mean"=mean(y),"median"=median(y),"variance"=var(y),"sd"=sd(y))
z


## IFELSE ##

ifelse(test,yes,no)

ifelse(x_scaled > 0, print(x_scaled), print("Negative"))
ifelse(x_scaled > 0, print(x_scaled), print(x_scaled^2))
as.vector(ifelse(x_scaled > 0, print(x_scaled), print(sqrt(x_scaled^2))))


## USER-WRITTEN FUNCTIONS ##

myfunction <- function(arg1, arg2, ... ){
statements
return(object)
}


square.it <- function(x) {
    square <- x * x
    return(square)
}
square.it(5)


# function example - get measures of central tendency
# and spread for a numeric vector x. The user has a
# choice of measures and whether the results are printed.
mysummary <- function(x,npar=TRUE,print=TRUE) {
  if (!npar) {
    center <- mean(x); spread <- sd(x) 
  } else {
    center <- median(x); spread <- mad(x) 
  }
  if (print & !npar) {
    cat("Mean=", center, "\n", "SD=", spread, "\n")
  } else if (print & npar) {
    cat("Median=", center, "\n", "MAD=", spread, "\n")
  }
  result <- list(center=center,spread=spread)
  return(result)
}

# invoking the function 
set.seed(1234)
x <- rpois(500, 4) 
y <- mysummary(x)
Median= 4
MAD= 1.4826 
# y$center is the median (4) 
# y$spread is the median absolute deviation (1.4826)

y <- mysummary(x, npar=FALSE, print=FALSE)
# no output 
# y$center is the mean (4.052)
# y$spread is the standard deviation (2.01927)

mysummary(x_scaled)
mysummary(x_scaled, npar=FALSE)
mysummary(x_scaled, print=FALSE)
mysummary(x_scaled, npar=FALSE, print=FALSE)
mysummary(x_scaled, npar=FALSE, print=FALSE)[["center"]]
mysummary(x_scaled, npar=FALSE, print=FALSE)[["spread"]]


## SORTING DATA ##

# sorting examples using the mtcars dataset
attach(mtcars)

# sort by mpg
mtcars[order(mpg),]
mtcars[order(-mpg),]
newdata <- mtcars[order(mpg),] 

# sort by mpg and cyl
newdata <- mtcars[order(mpg, cyl),]

#sort by mpg (ascending) and cyl (descending)
newdata <- mtcars[order(mpg, -cyl),]



## MERGING DATA ##

## ADDING COLUMNS ##

# merge two data frames by ID
total <- merge(dataframeA, dataframeB, by="ID")

# merge two data frames by ID and Country
total <- merge(dataframeA, dataframeB, by=c("ID","Country"))


d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
df1 <- data.frame(d,e,f)
names(df1) <- c("ID","CoatColor","Passed")

d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
df2 <- data.frame(d,e,f)
names(df2) <- c("ID","CarColor","Lied")

df1df2 <- merge(df1, df2, by="ID")
df1df2



## ADDING ROWS ##

# The two data frames must have the same variables, but they do not have to be in the same order
# If data frameA has variables that data frameB does not, then either:
# Delete the extra variables in data frameA or
# Create the additional variables in data frameB and set them to NA (missing) before joining them with rbind
total <- rbind(dataframeA, dataframeB)

names(df2) <- c("ID","CoatColor","Passed")
df1df2 <- rbind(df1, df2)
df1df2



## AGGREGATING DATA ##

# aggregate data frame mtcars by cyl and vs, returning means
# for numeric variables
attach(mtcars)
aggdata <- aggregate(mtcars, by=list(cyl,gear), FUN=mean, na.rm=TRUE)
print(aggdata)



## RESHAPING DATA ##

mydata
id	time	x1	x2
1	1	5	6
1	2	3	5
2	1	6	1
2	2	2	4

id <- c(1,1,2,2)
time <- c(1,2,1,2)
x1 <- c(5,3,6,2)
x2 <- c(6,5,1,4)
df3 <- data.frame(id, time, x1, x2)

# example of melt function 
install.packages("reshape")
library(reshape)
df3_melted <- melt(df3, id=c("id","time"), na.rm=TRUE)

df3_melted
id	time	variable	value
1	1	x1	5
1	2	x1	3
2	1	x1	6
2	2	x1	2
1	1	x2	6
1	2	x2	5
2	1	x2	1
2	2	x2	4

# cast the melted data
# cast(data, formula, function) 
subjmeans <- cast(df3_melted, id~variable, mean)
timemeans <- cast(df3_melted, time~variable, mean)

subjmeans
id	x1	x2
1	4	5.5
2	4	2.5

timemeans
time	x1	x2
1	5.5	3.5
2	2.5	4.5



## SUBSETTING DATA ##

## SELECTING (KEEPING) VARIABLES

# select variables time, x1
myvars <- c("time", "x1")
newdata <- df3[myvars]
newdata

# another method
myvars <- paste("v", 1:3, sep="")
newdata <- mydata[myvars]

# select 1st and 5th thru 10th variables
newdata <- mydata[c(1,5:10)]


## EXCLUDING (DROPPING) VARIABLES

# exclude variables time, x1
myvars <- names(df3) %in% c("time", "x1") 
newdata <- df3[!myvars]
newdata

# exclude 3rd and 5th variable 
newdata <- mydata[c(-3,-5)]

# delete variables v3 and v5
mydata$v3 <- mydata$v5 <- NULL


## SELECTING OBSERVATIONS ##

# first 5 observations
newdata <- mtcars[1:5,]
newdata

# based on variable values
newdata <- mtcars[ which(mtcars$cyl==6 & mtcars$mpg > 20), ]
newdata

# or
attach(newdata)
newdata <- mydata[ which(gender=='F' & age > 65),]
detach(newdata)


## SUBSET FUNCTION ##

# using subset function 
newdata <- subset(airquality, Ozone >= 100 | Solar.R < 10, select=c(Ozone, Solar.R, Wind, Temp, Month, Day))
newdata

# using subset function (part 2)
newdata <- subset(mydata, sex=="m" & age > 25,
select=weight:income)

newdata <- subset(airquality, Ozone >= 100 | Solar.R < 10, select=c(Ozone, Solar.R, Wind, Temp, Month, Day))
newdata


## RANDOM SAMPLE ##

# take a random sample of size 10 from a dataset airquality
# sample without replacement
mysample <- airquality[sample(1:nrow(airquality), 10, replace=FALSE),]
mysample



## DATA TYPE CONVERSION ##

Use is.foo to test for data type foo. Returns TRUE or FALSE
Use as.foo to explicitly convert it.

is.numeric(), is.character(), is.vector(), is.matrix(), is.data.frame()
as.numeric(), as.character(), as.vector(), as.matrix(), as.data.frame()



## WITH ##

# with(data, expression)
# example applying a t-test to a data frame mydata 
with(mydata, t.test(y ~ group))

with(df3, t.test(x2 ~ id))

with(mtcars[which(mtcars$cyl %in% c(4,8)), ], t.test(mpg ~ cyl))



## BY ##

# by(data, factorlist, function)
# example obtain variable means separately for
# each level of byvar in data frame mydata 
by(mydata, mydatat$byvar, function(x) mean(x))

by(mtcars, mtcars$cyl, function(x) mean(mtcars$mpg))

by(mtcars, mtcars$cyl, function(x) lm(mpg ~ wt, data=x))

output <- with(mtcars, by(mtcars, mtcars$cyl, function(x) lm(mpg ~ wt, data=x)))
sapply(output, coef)



## DESCRIPTIVE STATISTICS ##

# get means for variables in data frame mydata
# excluding missing values 
sapply(mydata, mean, na.rm=TRUE)
# possible functions used in sapply include mean, sd, var, min, max, median, range, and quantile

sapply(mtcars, mean, na.rm=TRUE)
sapply(airquality, range, na.rm=TRUE)


# mean,median,25th and 75th quartiles,min,max
summary(mtcars)
summary(airquality)


# Tukey min,lower-hinge, median,upper-hinge,max
fivenum(mtcars$mpg)


install.packages("Hmisc")
library(Hmisc)
# n, nmiss, unique, mean, 5,10,25,50,75,90,95th percentiles 
# 5 lowest and 5 highest scores
describe(airquality) 


## BY GROUP ##

# psych
install.packages("psych")
library(psych)
describe.by(mydata, group,...)

# doBy
install.packages("doBy")
library(doBy)
summaryBy(mpg + wt ~ cyl + vs, data = mtcars, 
 	FUN = function(x) { c(m = mean(x), s = sd(x)) } )
# produces mpg.m wt.m mpg.s wt.s for each 
# combination of the levels of cyl and vs



## FREQUENCIES AND CROSSTABS ##

# 2-Way Frequency Table 
attach(mtcars)
mytable <- table(cyl,gear) # A will be rows, B will be columns 
mytable # print table 

margin.table(mytable, 1) # A frequencies (summed over B) 
margin.table(mytable, 2) # B frequencies (summed over A)

prop.table(mytable) # cell percentages
prop.table(mytable, 1) # row percentages 
prop.table(mytable, 2) # column percentages



## CORRELATIONS ##

A simplified format is cor(x, use=, method= ) where

Option	Description
x	Matrix or data frame

use	Specifies the handling of missing data. Options are all.obs (assumes no missing data - missing data will produce an error), complete.obs (listwise deletion), and pairwise.complete.obs (pairwise deletion)

method	Specifies the type of correlation. Options are pearson, spearman or kendall.

# Correlations/covariances among numeric variables in 
# data frame mtcars. Use listwise deletion of missing data. 
cor(mtcars, use="complete.obs", method="kendall") 
cov(mtcars, use="complete.obs")


# Correlations with significance levels
# Input must be a matrix; Uses pairwise deletion
install.packages("Hmisc")
library(Hmisc)
rcorr(x, type="pearson") # type can be pearson or spearman

#mtcars is a data frame
rcorr(as.matrix(mtcars))

cor(mtcars, use="complete.obs", method="pearson") 



## GRAPHING AND PLOTTING ##

## CREATING A GRAPH ##
attach(mtcars)
plot(wt, mpg) 
abline(lm(mpg~wt))
title("Regression of MPG on Weight")


## SAVING A GRAPH ##

Function	Output to
pdf("mygraph.pdf")	# pdf file
win.metafile("mygraph.wmf")	windows metafile
png("mygraph.png")	png file
jpeg("mygraph.jpg")	jpeg file
bmp("mygraph.bmp")	bmp file
postscript("mygraph.ps")	postscript file
dev.off()

# For example
# High Density Scatterplot with Color Transparency 
pdf("c:/scatterplot.pdf") 
x <- rnorm(1000)
y <- rnorm(1000) 
plot(x,y, main="PDF Scatterplot Example", col=rgb(0,100,0,50,maxColorValue=255), pch=16)
dev.off()


## HISTOGRAMS ##

# Simple Histogram
hist(mtcars$mpg)

# Colored Histogram with Different Number of Bins
hist(mtcars$mpg, breaks=12, col="red")

# Add a Normal Curve (Thanks to Peter Dalgaard)
x <- mtcars$mpg 
h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon", main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)


## DENSITY PLOTS ##

# Kernel Density Plot
d <- density(mtcars$mpg) # returns the density data 
plot(d) # plots the results

# Filled Density Plot
d <- density(mtcars$mpg)
plot(d, main="Kernel Density of Miles Per Gallon")
polygon(d, col="red", border="blue")


## DOT PLOTS ##

# Simple Dotplot
dotchart(mtcars$mpg,labels=row.names(mtcars),cex=.7,
  	main="Gas Milage for Car Models", 
   xlab="Miles Per Gallon")

# Dotplot: Grouped Sorted and Colored
# Sort by mpg, group and color by cylinder 
x <- mtcars[order(mtcars$mpg),] # sort by mpg
x$cyl <- factor(x$cyl) # it must be a factor
x$color[x$cyl==4] <- "red"
x$color[x$cyl==6] <- "blue"
x$color[x$cyl==8] <- "darkgreen"	
dotchart(x$mpg,labels=row.names(x),cex=.7,groups= x$cyl,
  	main="Gas Milage for Car Models\ngrouped by cylinder",
   xlab="Miles Per Gallon", gcolor="black", color=x$color)


## BAR PLOTS ##

# Simple Bar Plot 
counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution", xlab="Number of Gears")

# Simple Horizontal Bar Plot with Added Labels 
counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution", horiz=TRUE, names.arg=c("3 Gears", "4 Gears", "5 Gears"))

# Stacked Bar Plot with Colors and Legend
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
  xlab="Number of Gears", col=c("darkblue","red"),
 	legend = rownames(counts))

# Grouped Bar Plot
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
  xlab="Number of Gears", col=c("darkblue","red"),
 	legend = rownames(counts), beside=TRUE)

# NOTE: Bar plots need not be based on counts or frequencies. You can create bar plots that represent means, medians, standard deviations, etc. Use the aggregate function and pass the results to the barplot function


## LINE CHARTS ##

type	description
p	points
l	lines
o	overplotted points and lines
b,c	points (empty if "c") joined by lines
s,S	stair steps
h	histogram-like vertical lines
n	does not produce any points or lines

x <- c(1:5); y <- x # create some data 
par(pch=22, col="red") # plotting symbol and color 
par(mfrow=c(2,4)) # all plots on one page 
opts = c("p","l","o","b","c","s","S","h") 
for(i in 1:length(opts)){ 
  heading = paste("type=",opts[i]) 
  plot(x, y, type="n", main=heading) 
  lines(x, y, type=opts[i]) 
}

x <- c(1:5); y <- x # create some data
par(pch=22, col="blue") # plotting symbol and color 
par(mfrow=c(2,4)) # all plots on one page 
opts = c("p","l","o","b","c","s","S","h") 
for(i in 1:length(opts)){ 
  heading = paste("type=",opts[i]) 
  plot(x, y, main=heading) 
  lines(x, y, type=opts[i]) 
}


## BOX PLOTS ##

# Boxplot of MPG by Car Cylinders 
boxplot(mpg~cyl,data=mtcars, main="Car Milage Data", xlab="Number of Cylinders", ylab="Miles Per Gallon")

# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation 
boxplot(len~supp*dose, data=ToothGrowth, notch=TRUE, 
  col=(c("gold","darkgreen")),
  main="Tooth Growth", xlab="Suppliment and Dose")

# NOTE: In the notched boxplot, if two boxes' notches do not overlap this is ‘strong evidence’ their medians differ


## SCATTER PLOTS ##

# Simple Scatterplot
attach(mtcars)
plot(wt, mpg, main="Scatterplot Example", xlab="Car Weight ", ylab="Miles Per Gallon ", pch=19)
# Add fit lines
abline(lm(mpg~wt), col="red") # regression line (y~x) 
lines(lowess(wt,mpg), col="blue") # lowess line (x,y)

# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
scatterplot(mpg ~ wt | cyl, data=mtcars, 
  	xlab="Weight of Car", ylab="Miles Per Gallon", 
   main="Enhanced Scatter Plot", 
   labels=row.names(mtcars))

# Basic Scatterplot Matrix
pairs(~mpg+disp+drat+wt,data=mtcars, 
   main="Simple Scatterplot Matrix")

# Scatterplot Matrices from the car Package
scatterplot.matrix(~mpg+disp+drat+wt|cyl, data=mtcars,
  	main="Three Cylinder Options")

# High Density Scatterplot with Binning
install.packages("hexbin")
library(hexbin)
x <- rnorm(1000)
y <- rnorm(1000)
bin<-hexbin(x, y, xbins=50) 
plot(bin, main="Hexagonal Binning")

# High Density Scatterplot with Color Transparency 
pdf("c:/scatterplot.pdf") 
x <- rnorm(1000)
y <- rnorm(1000) 
plot(x,y, main="PDF Scatterplot Example", col=rgb(0,100,0,50,maxColorValue=255), pch=16)
dev.off()

# 3D Scatterplot
install.packages("scatterplot3d")
library(scatterplot3d)
attach(mtcars)
scatterplot3d(wt,disp,mpg, main="3D Scatterplot")

# 3D Scatterplot with Coloring and Vertical Drop Lines
library(scatterplot3d) 
attach(mtcars) 
scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE,
  type="h", main="3D Scatterplot")

# 3D Scatterplot with Coloring and Vertical Lines
# and Regression Plane 
library(scatterplot3d) 
attach(mtcars) 
s3d <-scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE,
  type="h", main="3D Scatterplot")
fit <- lm(mpg ~ wt+disp) 
s3d$plane3d(fit)
