# authors   Juris Freimanis and Romans Pezenkovs

import matplotlib.pyplot
import math
import numpy as np

image_x_dim = 200
image_y_dim = 200
cell_size = 10 # size of cell for polarization vectors (if 10 then 10x10=100 )


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


def summ_several_cells(big_2d_list):
    short_2d_list = [[0.0 for x in range(int(image_x_dim/cell_size))] for y in range(int(image_y_dim/cell_size))]
    for y in range(int(image_y_dim/cell_size) ):
        for x in range(int(image_x_dim/cell_size) ):
            for y_2 in range(cell_size):
                for x_2 in range(cell_size):
                    short_2d_list[y][x] += big_2d_list[y*cell_size+y_2][x*cell_size+x_2]
    return short_2d_list


image_data_2d_list = init_2d_list('data_gnu_I_2.txt')
image_data_2d_list = np.array(image_data_2d_list)*1.874823e-31
data_for_vectors = summ_several_cells(image_data_2d_list)

min_val = -0.2*3e-30
max_val = 3e-30

for i in range(image_y_dim):
    for j in range(image_x_dim):
        if(image_data_2d_list[i][j] == 0.0):
            image_data_2d_list[i][j] = min_val

for i in range(int(image_y_dim/cell_size) ):
    for j in range(int(image_x_dim/cell_size) ):
        if(data_for_vectors[i][j] == 0.0):
            data_for_vectors[i][j] = min_val

data_Q2_list = init_2d_list('data_gnu_Q_2.txt')
data_U2_list = init_2d_list('data_gnu_U_2.txt')
data_Q2_list = np.array(data_Q2_list)*1.874823e-31
data_U2_list = np.array(data_U2_list)*1.874823e-31

data_Q2 = summ_several_cells(data_Q2_list)
data_U2 = summ_several_cells(data_U2_list)

min_pol_degree = 1 # minimum value of polarization degree
max_pol_degree = 0 # maximum value of polarization degree

for y in range( image_y_dim ):
    for x in range( image_x_dim ):
        if(image_data_2d_list[y][x] != min_val):
            pol_v_len = math.sqrt(data_Q2_list[y][x]**2 + data_U2_list[y][x]**2) / image_data_2d_list[y][x]
            if(pol_v_len < min_pol_degree):
                min_pol_degree = pol_v_len
            if(pol_v_len > max_pol_degree):
                max_pol_degree = pol_v_len

print('Minimum polarization degree = ', min_pol_degree)
print('Maximum polarization degree = ', max_pol_degree)
cmap_val = 'viridis'
#cmap_val = 'cividis'
#cmap_val = 'inferno'
#cmap_val = 'summer'

matplotlib.pyplot.imshow(image_data_2d_list, cmap = cmap_val, vmin = min_val, vmax = max_val, origin = 'lower' )

my_colorbar = matplotlib.pyplot.colorbar()
my_colorbar.set_label(r'$H^{pix}_{\nu}$' + '  [W / Hz]')
ax = my_colorbar.ax
text = ax.yaxis.label
font = matplotlib.font_manager.FontProperties(size=14)
text.set_font_properties(font)

matplotlib.pyplot.xlabel('$\\lambda$' + ' = 550nm')

matplotlib.pyplot.title('Integral_equ, opt_rad=0.1, albedo=0.99')

max_pol_degree = 0 # maximum value of polarization degree

for y in range(int(image_y_dim/cell_size) ):
    for x in range(int(image_x_dim/cell_size) ):
        if(data_for_vectors[y][x] != min_val):
            pol_v_len = math.sqrt(data_Q2[y][x]*data_Q2[y][x]+data_U2[y][x]*data_U2[y][x]) / data_for_vectors[y][x]                
            if(pol_v_len > max_pol_degree):
                max_pol_degree = pol_v_len
            angle_polar = math.atan2( -data_U2[y][x], data_Q2[y][x]) / 2.0
            matplotlib.pyplot.arrow(x*cell_size + 5, y*cell_size + 5, 5*pol_v_len*math.cos(angle_polar), 5*pol_v_len*math.sin(angle_polar), width = 0.001, edgecolor='white', facecolor='white')
            matplotlib.pyplot.arrow(x*cell_size + 5, y*cell_size + 5, -5*pol_v_len*math.cos(angle_polar), -5*pol_v_len*math.sin(angle_polar), width = 0.001, edgecolor='white', facecolor='white')
matplotlib.pyplot.arrow(185, 185, 5*max_pol_degree, 0, width = 0.001, edgecolor='white', facecolor='white')
matplotlib.pyplot.arrow(185, 185, -5*max_pol_degree, 0, width = 0.001, edgecolor='white', facecolor='white')
str_pol_deg = 'p=' + str(round(max_pol_degree*100, 1)) + '%'
matplotlib.pyplot.text(165, 190, str_pol_deg , {'color': 'w', 'fontsize': 10})
print('10x10 cells Maximum polarization degree = ', max_pol_degree)

matplotlib.pyplot.show()

