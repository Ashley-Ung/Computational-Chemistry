"""
txtFileToExcel.py
Ashley Ung
Spring 2023

This Python program extracts the numbers (in the row called "averages" that are space delimitated) from two files called filename_hydrophobic_interface_comp.txt and filename_aqueous_interface_comp.txt where "filename" is some file that was read in for analysis and will change for each time the code runs and places the averages into an Excel file called averages.xlsx.
Must install the xlwt library to create Excel files and pandas to read the text files and extract the row called "averages". pip3 install xlwt, pip3 install pandas
This code assumes that the text files are formatted in a way that can be loaded by pandas, and that the column name for the averages row is "Averages" and has to be in the second column. used the "iloc" property of the Dataframe to extract the second column, which is the column containing the averages.
"""
import xlwt
import pandas as pd

# Read the filename from user input
filename = input ("Enter the filename: ")

# Read data from filename_hydrophobic_interface_comp.txt and filename_aqueous_interface_comp.txt
data1 = pd.read_csv (filename + '_hydrophobic_interface_comp.txt', delimiter=' ')
data2 = pd.read_csv (filename + '_aqueous_interface_comp.txt', delimiter=' ')

# Extract averages from the second column
avg1 = data1.iloc[:, 1]
avg2 = data2.iloc[:, 1]

# Create workbook and worksheet
workbook = xlwt.Workbook ()
worksheet = workbook.add_sheet ('Averages')

# Write data to worksheet
for i in range (len (avg1)):
    worksheet.write (i, 0, avg1[i])
    worksheet.write (i, 1, avg2[i])
    
# Save workbook
workbook.save ('averages.xlsx')
