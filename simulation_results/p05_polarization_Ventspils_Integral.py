# authors   Juris Freimanis and Romans Pezenkovs

import matplotlib.pyplot
import math

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


data_I_1 = init_2d_list('data_gnu_I.txt')
data_I_2 = init_2d_list('data_gnu_I_2.txt')

data_Q_1 = init_2d_list('data_gnu_Q.txt')
data_U_1 = init_2d_list('data_gnu_U.txt')

data_Q_2 = init_2d_list('data_gnu_Q_2.txt')
data_U_2 = init_2d_list('data_gnu_U_2.txt')

compare_2d_list = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)]


for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(data_I_2[y][x] != 0 and data_I_1[y][x] != 0 and math.sqrt(data_Q_2[y][x]**2 + data_U_2[y][x]**2) != 0 ):
            compare_2d_list[y][x] = (math.sqrt(data_Q_1[y][x]**2 + data_U_1[y][x]**2) / data_I_1[y][x] ) / \
                                    (math.sqrt(data_Q_2[y][x]**2 + data_U_2[y][x]**2) / data_I_2[y][x] )

number_of_elements = 0
summ_of_all_elements = 0

for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        if(compare_2d_list[y][x] > 0):
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
print('standard_deviation = ', standard_deviation)
print('relative_error = ', relative_error)

            
matplotlib.pyplot.imshow(compare_2d_list, cmap = 'inferno', vmin = 0, vmax = 2,  origin = 'lower' )
my_colorbar = matplotlib.pyplot.colorbar()

my_colorbar.set_label('Ventspils / Integral_equ')
ax = my_colorbar.ax
text = ax.yaxis.label
#font = matplotlib.font_manager.FontProperties(family='times new roman', style='italic', size=16)
font = matplotlib.font_manager.FontProperties(size=12)
text.set_font_properties(font)
matplotlib.pyplot.xlabel('$\\lambda$' + ' = 550nm')

matplotlib.pyplot.title('Polarization deg., opt_way=0.1, albedo=0.99')
str_average_val = '1.018' + u'\u00B1' + '0.013'
matplotlib.pyplot.text(3, 6, str_average_val , {'color': 'black', 'fontsize': 12}, bbox=dict(facecolor='w', alpha=1))
matplotlib.pyplot.show()

