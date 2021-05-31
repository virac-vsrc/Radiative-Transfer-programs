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
data_Q2_list = [[0.0 for x in range(image_x_dim)] for y in range(image_y_dim)]
data_U2_list = [[0.0 for x in range(image_x_dim)] for y in range(image_y_dim)]

y_radmc = 0

for y in range(image_y_dim):
    for x in range(image_x_dim):
        if(image_out[y_radmc][0] > 1e-40):
            image_data_2d_list[y][x] = image_out[y_radmc][0]
            data_Q2_list[y][x] = image_out[y_radmc][1]
            data_U2_list[y][x] = -image_out[y_radmc][2] # U=-U RADMC and M.I.Mishchenko 
        y_radmc = int(y_radmc+1)


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


data_Q2 = summ_several_cells(data_Q2_list)
data_U2 = summ_several_cells(data_U2_list)

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

matplotlib.pyplot.title('RADMC3D, opt_rad=0.1, albedo=0.99')

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
matplotlib.pyplot.show()
