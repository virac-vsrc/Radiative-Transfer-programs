# authors   Juris Freimanis and Romans Pezenkovs

import matplotlib.pyplot
import math
import numpy as np

image_x_dim = 200
image_y_dim = 200


def init_2d_list(file_name):
    data_2d_list = [[0.0 for x in range(image_x_dim)] for y in range(image_y_dim)]
    # to access to variable of 2d list data_2d_list[y][x]
    data_file = open(file_name)
    row_numb = 0
    while True:
        line = data_file.readline()
        if(len(line) == 0):
            break
        number_list = [float(k) for k in line[:-2].split(' ')] # line[:-2] except two last symbols, split numbers when symbol is ' '
        data_2d_list[row_numb] = number_list
        row_numb += 1
    data_file.close()
    return data_2d_list


image_data_2d_list_1 = init_2d_list('data_gnu_I.txt')
image_data_2d_list_2 = init_2d_list('data_gnu_I_2.txt')
image_data_2d_list_1 = np.array(image_data_2d_list_1)*4.997579e-01
image_data_2d_list_2 = np.array(image_data_2d_list_2)*1.874823e-31
compare_2d_list = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)]

number_of_elements = 0
summ_of_all_elements = 0

for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(image_data_2d_list_2[y][x] != 0):
            compare_2d_list[y][x] = image_data_2d_list_1[y][x] / image_data_2d_list_2[y][x]
            number_of_elements += 1
            summ_of_all_elements += compare_2d_list[y][x]

average = summ_of_all_elements / number_of_elements

summ_standard_deviation = 0

for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(image_data_2d_list_2[y][x] != 0):
            summ_standard_deviation += (compare_2d_list[y][x] - average)**2

standard_deviation = math.sqrt( summ_standard_deviation / (number_of_elements - 1) )

relative_error = standard_deviation / average

print('average = ', average)
print('standard_deviation = ', standard_deviation)
print('relative_error = ', relative_error)
min_v = 0
max_v = 2

matplotlib.pyplot.imshow(compare_2d_list, cmap = 'inferno', vmin = min_v ,vmax = max_v, origin = 'lower' )

my_colorbar = matplotlib.pyplot.colorbar()
my_colorbar.set_label('Ventspils / Integral_equ')
ax = my_colorbar.ax
text = ax.yaxis.label
#font = matplotlib.font_manager.FontProperties(family='times new roman', style='italic', size=16)
font = matplotlib.font_manager.FontProperties(size=12)
text.set_font_properties(font)
matplotlib.pyplot.xlabel('$\\lambda$' + ' = 550nm')

matplotlib.pyplot.title('Intensity, opt_rad=0.1, albedo=0.99')
str_average_val = '0.965' + u'\u00B1' + '0.130'
matplotlib.pyplot.text(3, 6, str_average_val , {'color': 'black', 'fontsize': 12}, bbox=dict(facecolor='w', alpha=1))
matplotlib.pyplot.show()

