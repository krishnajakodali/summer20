import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlsxwriter
plt.rcdefaults()
#reading the table
loc=("/home/krishnajakodali/probstat/stat/tables/exer40/cont40.xlsx")
wb=xlrd.open_workbook(loc)
sheet=wb.sheet_by_index(0)
#grouped frequency distribution table specifications
class_width=9
start_value=117.5
end_value=180.5
#number of classes
total_classes=int((end_value-start_value)/class_width)
#combining all data
data=[]
for i in range(total_classes):
	for j in range(int(sheet.cell_value(i+1,1))):
		data.append(start_value+i*class_width+1)
#creating histogram
#giving bin_edges
bin_edges=[]
for m in range(total_classes+1):
	bin_edges.append(0)
for x in range(total_classes+1):
	bin_edges[x]=start_value+x*class_width
#plotting the histogram
plt.hist(data,bins=bin_edges)
plt.xlabel("Length of leaves in mm")
plt.ylabel("number of leaves")
plt.savefig('./pyfigs/exer40.eps')


plt.show()

