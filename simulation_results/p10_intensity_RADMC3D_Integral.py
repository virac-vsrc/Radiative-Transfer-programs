# authors   Juris Freimanis and Romans Pezenkovs

import matplotlib.pyplot
import math
import numpy as np

image_x_dim = 200
image_y_dim = 200
cell_size = 10 # size of cell for polarization vectors (if 10 then 10x10=100 )

radmc_x_dim = 4 # 4 Stokes parameters
radmc_y_dim = int(image_x_dim*image_y_dim)

def read_image_out(file_name):
    data_2d_list = [[0.0 for x in range(radmc_x_dim)] for y in range(radmc_y_dim)]
    data_file = open(file_name)
    row_numb = 0
    while True:
        line = data_file.readline()
        if(len(line) == 0):
            break
        nums = line.split()
        nums = list(map(float, nums))
        data_2d_list[row_numb] = nums
        row_numb += 1
    data_file.close()
    return data_2d_list


def init_2d_list(file_name):
    data_2d_list = [[0.0 for x in range(image_x_dim)] for y in range(image_y_dim)]
    # to access to variable of 2d list data_2d_list[y][x]
    data_file = open(file_name)
    row_numb = 0
    while True:
        line = data_file.readline()
        if(len(line) == 0):
            break
        nums = line.split()
        nums = list(map(float, nums))
        data_2d_list[row_numb] = nums
        row_numb += 1
    data_file.close()
    return data_2d_list


def summ_several_cells(big_2d_list):
    short_2d_list = [[0.0 for x in range(int(image_x_dim/cell_size))] for y in range(int(image_y_dim/cell_size))]
    for y in range(int(image_y_dim/cell_size) ):
        for x in range(int(image_x_dim/cell_size) ):
            for y_2 in range(cell_size):
                for x_2 in range(cell_size):
                    short_2d_list[y][x] += big_2d_list[y*cell_size+y_2][x*cell_size+x_2]
    return short_2d_list


image_out = read_image_out("image.out")
image_out = np.array(image_out)*1.436060e-19

image_data_2d_list = [[0.0 for x in range(image_x_dim)] for y in range(image_y_dim)]

y_radmc = 0

for y in range(image_y_dim):
    for x in range(image_x_dim):
         image_data_2d_list[y][x] = image_out[y_radmc][0]
         y_radmc = int(y_radmc+1)


data_I_1 = init_2d_list('data_gnu_I_2.txt')
data_I_1 = np.array(data_I_1)*1.874823e-31
compare_2d_list = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)]

number_of_elements = 0
summ_of_all_elements = 0

for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(image_data_2d_list[y][x] > 1e-40 and data_I_1[y][x]> 1e-40):
            if(image_data_2d_list[y][x] / data_I_1[y][x] < 10):
                compare_2d_list[y][x] = image_data_2d_list[y][x] / data_I_1[y][x] 
                number_of_elements += 1
                summ_of_all_elements += compare_2d_list[y][x]

average = summ_of_all_elements / number_of_elements

summ_standard_deviation = 0

for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(compare_2d_list[y][x] > 0):
            summ_standard_deviation += ( compare_2d_list[y][x] - average )**2

standard_deviation = math.sqrt( summ_standard_deviation / (number_of_elements - 1) )            

relative_error = standard_deviation / average

print('average = ', average)
print('{:.4e}'.format(average) )
print('standard_deviation = ', standard_deviation)
print('{:.4e}'.format(standard_deviation) )
print('relative_error = ', relative_error)
print('{:.4e}'.format(relative_error) )

matplotlib.pyplot.imshow(compare_2d_list, cmap = 'inferno', vmin = 0, vmax = 2, origin = 'lower' )

my_colorbar = matplotlib.pyplot.colorbar()
my_colorbar.set_label('RADMC3D / Integral_equ')
ax = my_colorbar.ax
text = ax.yaxis.label

font = matplotlib.font_manager.FontProperties(size=12)
text.set_font_properties(font)
matplotlib.pyplot.xlabel('$\\lambda$' + ' = 550nm')

matplotlib.pyplot.title('Intensity, opt_rad=0.1, alb.=0.99')
str_average_val = '1.036' + u'\u00B1' + '0.063'
matplotlib.pyplot.text(3, 6, str_average_val , {'color': 'black', 'fontsize': 12}, bbox=dict(facecolor='w', alpha=1))

matplotlib.pyplot.show()
