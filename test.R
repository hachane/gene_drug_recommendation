library(readxl)
library(ComplexHeatmap)
library(dplyr)
# Read the Excel file
table1c <- read_xlsx("allele frequency.xlsx")

# Extract the numeric data (excluding the 'Gene' column)
numeric_data <- as.matrix(table1c[, -1])

rownames(numeric_data) <- numeric_data[, 1]

# Remove the first column from the numeric data
numeric_data <- numeric_data[, -1]
index_array <- rownames(numeric_data)

# Extract column names
column_array <- colnames(numeric_data)
#convert the datatype of "mat" from "char" to "numeric" number
numeric_data = matrix(as.numeric(numeric_data),    ncol = ncol(numeric_data))
rownames(numeric_data) <- index_array
colnames(numeric_data) <- column_array
Heatmap(numeric_data)