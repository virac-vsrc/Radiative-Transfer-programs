# authors   Juris Freimanis and Romans Pezenkovs

import math
from scipy import special
import numpy as np

optical_radius = 0.1 # optical radius of the cloud only
                   # full optical radius is 'optical_radius + tau_star'

albedo = 0.4664 # albedo of dust cloud

tau = 0
tau_prime = 0
tau_star = 0.00256 * optical_radius # optical radius of the star
                                   # tau_star should be less than tau and tau_prime

t_max = 0 # parameter for formula (14)

full_opt_rad = optical_radius + tau_star # full optical radius of the cloud

mu = 0 # cosine of the angle between the radial direction and the direction of propagation
       # of radiation

Luminosity = 3.828e+26 # star's luminosity
                       # Luminosity of sun 3.828e+26 Watt

parsec = 3.0856775815 * 10**16     # one parsec in meters
astron_unit = 1.495978707 * 10**11 # one astronomical unit in meters

cloud_radius = 20 * astron_unit # radius of cloud in meters
dist_cloud_tel = 600 * parsec   # distance between cloud and telescope
dist_star_cent_tel = cloud_radius + dist_cloud_tel # distance between star center and telescope

a = full_opt_rad / cloud_radius # extinction coefficient, to convert geometrical radius to optical

width = 4001    # Point number for tau_prime, should be odd number for Simpsons method 
height = width  # Point number for tau, should be the same as width

image_x_dim = 200           # The number of pixels (in X direction) for image of light scattered by dust cloud
image_y_dim = image_x_dim   # The number of pixels (in Y direction) for image of light scattered by dust cloud
                            # if 200x200 then 40000 pixels in calculated image

K00_array = [[0 for x in range(width)] for y in range(height)] # kernel K00 of the equation (12)
K02_array = [[0 for x in range(width)] for y in range(height)] # kernel K02 of the equation (12)
K20_array = [[0 for x in range(width)] for y in range(height)] # kernel K20 of the equation (13)
K22_array = [[0 for x in range(width)] for y in range(height)] # kernel K22 of the equation (13)
B00_array = [0 for y in range(height)] # B00, the expression (20)
B20_array = [0 for y in range(height)] # B20, the expression (21)
f02_array = [0 for y in range(height)] # f02 is a function bounded at tau = tau_prime, used in equation (43)
f20_array = [0 for y in range(height)] # f20 is a function bounded at tau = tau_prime, used in equation (48)
f22_array = [0 for y in range(height)] # f22 is a function bounded at tau = tau_prime, used in equation (48)
I00_array = [0 for y in range(height)] # I00, the expression (40)
I02_array = [0 for y in range(height)] # I02, the expression (42)
I20_array = [0 for y in range(height)] # I20, the expression (45)
I22_array = [0 for y in range(height)] # I22, the expression (47)
B0_array = [0 for y in range(height)] # B0 is an unknown function of the equations (12),(13)
B2_array = [0 for y in range(height)] # B2 ia an unknown function of the equations (12),(13)
I_array = [0 for y in range(height)] # Stokes parameter I equation (2)
Q_array = [0 for y in range(height)] # Stokes parameter Q equation (3)
data_gnu_I_2 = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)] # two dimensional list of Stokes parameter I for image
data_gnu_Q_2 = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)] # two dimensional list of Stokes parameter Q for image
data_gnu_U_2 = [[0 for x in range(image_x_dim)] for y in range(image_y_dim)] # two dimensional list of Stokes parameter U for image

delta_1 = 0.5 / (height-1) * optical_radius
delta_2 = 0.5 / (height-1) * optical_radius

for tau_p in range(0, height, 1): # K00, K02, K20, K22
    for tau_prime_p in range(0, width, 1):
        tau = tau_p / (height-1) * optical_radius
        tau_prime = tau_prime_p / (width-1) * optical_radius
        tau += tau_star
        tau_prime += tau_star

        # expression (14)
        t_max = math.sqrt(tau**2 - tau_star**2) + math.sqrt(tau_prime**2 - tau_star**2)

        if(tau == tau_star and tau_prime == tau_star):
            K00_array[tau_p][tau_prime_p] = 0
            K02_array[tau_p][tau_prime_p] = 0
            K20_array[tau_p][tau_prime_p] = 0
            K22_array[tau_p][tau_prime_p] = 0
        else:
            # expression (15)
            K00_array[tau_p][tau_prime_p] = albedo*tau_prime/(2*tau) * \
                ( special.exp1(abs(tau-tau_prime)) - special.exp1(t_max))

            # expression (16)
            K02_array[tau_p][tau_prime_p] = 3*albedo/(16*tau*tau_prime)* \
                ( math.exp(-abs(tau-tau_prime))*( abs(tau-tau_prime) + 1 + \
                (tau+tau_prime)**2/2 * (1-abs(tau-tau_prime)) ) - \
                math.exp(-t_max)*(t_max+1+(tau**2-tau_prime**2)**2*(1-t_max)/(2*t_max**2)) + \
                ( 2/3*tau_prime**2 - 2*tau**2 + (tau**2-tau_prime**2)**2/2 ) * \
                ( special.exp1(abs(tau-tau_prime)) - special.exp1(t_max)))

            # expression (17)
            K20_array[tau_p][tau_prime_p] = 3*albedo*tau_prime/(32*tau**3)* \
                ( math.exp(-abs(tau-tau_prime))*( abs(tau-tau_prime) + 1 + \
                (tau+tau_prime)**2/2 * (1-abs(tau-tau_prime)) ) - \
                math.exp(-t_max)*(t_max+1+(tau**2-tau_prime**2)**2*(1-t_max)/(2*t_max**2)) + \
                ( 2/3*tau**2 - 2*tau_prime**2 + (tau**2-tau_prime**2)**2/2 ) * \
                ( special.exp1(abs(tau-tau_prime)) - special.exp1(t_max)))

            # expression (18)
            K22_array[tau_p][tau_prime_p] = albedo/(128*tau**3*tau_prime)* \
                ( math.exp(-abs(tau-tau_prime)) * ( 9*abs(tau-tau_prime)**3 + \
                27*(tau-tau_prime)**2 + 54*abs(tau-tau_prime) + 54 - \
                24*(tau**2+tau_prime**2)*(abs(tau-tau_prime)+1) + \
                12*(tau+tau_prime)**2*(tau**2 + tau_prime**2)*(abs(tau-tau_prime)-1) + \
                3*(tau+tau_prime)**4/8 * ( 6-2*abs(tau-tau_prime)+(tau-tau_prime)**2 - \
                abs(tau - tau_prime)**3 ) ) - math.exp(-t_max) * ( 9*t_max**3 + \
                27*t_max**2 + 54*t_max + 54 - 24*(tau**2+tau_prime**2)*(t_max+1) + \
                12*(tau**2-tau_prime**2)**2*(tau**2+tau_prime**2)*(t_max-1)/t_max**2 + \
                3*(tau**2-tau_prime**2)**4/(8*t_max**4)*(6-2*t_max+t_max**2-t_max**3) ) + \
                ( 30*tau**4 + 20*tau**2*tau_prime**2 + 30*tau_prime**4 - \
                12*(tau**2-tau_prime**2)**2*(tau**2+tau_prime**2) + \
                3*(tau**2-tau_prime**2)**4/8 ) * \
                ( special.exp1(abs(tau-tau_prime)) - special.exp1(t_max) ) )

print('Kernels of integral equations K00, K02, K20, K22 are calculated.')

for tau_p in range(0, height, 1): # B00, B20
    tau = tau_p / (height-1) * optical_radius
    tau += tau_star

    if(tau == tau_star):
        B00_array[tau_p] = albedo*Luminosity*a**2/(8*math.pi**2*tau_star**2)
        B20_array[tau_p] = 0
    else:
        # expression (20)
        B00_array[tau_p] = albedo*Luminosity*a**2/(16*math.pi**2*tau_star**2*tau) * \
            ( (tau+tau_star-1)*math.exp(-(tau-tau_star)) - \
            ( math.sqrt(tau**2 - tau_star**2)-1 ) * math.exp(-math.sqrt(tau**2-tau_star**2)) + \
            (tau**2 - tau_star**2) * ( special.exp1(math.sqrt(tau**2-tau_star**2)) - \
            special.exp1(tau-tau_star) ) )

        # expression (21)
        B20_array[tau_p] = albedo*Luminosity*a**2/(256*math.pi**2*tau_star**2*tau**3) * \
            ( math.exp(-math.sqrt(tau**2-tau_star**2)) * ( (tau**2-tau_star**2)**2/2 + \
            2*tau**2 - 6*tau_star**2 + 6 + math.sqrt(tau**2-tau_star**2) * \
            ( 4*tau_star**2 + 6 - (tau**2-tau_star**2)**2/2 ) ) + math.exp(-(tau-tau_star)) * \
            ( (tau**2 - tau_star**2)*(tau+tau_star)**2*(tau-tau_star-1)/2 + \
            +2*tau_star*(tau**2-tau_star**2)-2*tau**2+6*tau*tau_star-6*(tau-tau_star)-6 ) + \
            (tau**2 - tau_star**2) * ( (tau**2-tau_star**2)**2/2 - tau**2 - 3*tau_star**2 ) * \
            ( special.exp1(math.sqrt(tau**2-tau_star**2)) - special.exp1(tau-tau_star) ) )



for tau_p in range(0, height, 1): # f02, f20, f22
    tau = tau_p / (height-1) * optical_radius
    tau += tau_star
    tau_prime = tau
    t_max = math.sqrt(tau**2 - tau_star**2) + math.sqrt(tau_prime**2 - tau_star**2)

    if(tau == tau_star):
        f02_array[tau_p] = 0
        f20_array[tau_p] = 0
        f22_array[tau_p] = 0
    else:

        # formula (41)
        f02_array[tau_p] = 3*albedo/(16*tau*tau_prime)* \
            ( math.exp(-abs(tau-tau_prime))*( abs(tau-tau_prime) + 1 + \
            (tau+tau_prime)**2/2 * (1-abs(tau-tau_prime)) ) - \
            math.exp(-t_max)*(t_max+1+(tau**2-tau_prime**2)**2*(1-t_max)/(2*t_max**2)) + \
            ( 2/3*tau_prime**2 - 2*tau**2 + (tau**2-tau_prime**2)**2/2 ) * \
            ( -special.exp1(t_max)))

        # formula (44)
        f20_array[tau_p] = 3*albedo*tau_prime/(32*tau**3)* \
            ( math.exp(-abs(tau-tau_prime))*( abs(tau-tau_prime) + 1 + \
            (tau+tau_prime)**2/2 * (1-abs(tau-tau_prime)) ) - \
            math.exp(-t_max)*(t_max+1+(tau**2-tau_prime**2)**2*(1-t_max)/(2*t_max**2)) + \
            ( 2/3*tau**2 - 2*tau_prime**2 + (tau**2-tau_prime**2)**2/2 ) * \
            ( -special.exp1(t_max)))

        # formula (46)
        f22_array[tau_p] = albedo/(128*tau**3*tau_prime)* \
            ( math.exp(-abs(tau-tau_prime)) * ( 9*abs(tau-tau_prime)**3 + \
            27*(tau-tau_prime)**2 + 54*abs(tau-tau_prime) + 54 - \
            24*(tau**2+tau_prime**2)*(abs(tau-tau_prime)+1) + \
            12*(tau+tau_prime)**2*(tau**2 + tau_prime**2)*(abs(tau-tau_prime)-1) + \
            3*(tau+tau_prime)**4/8 * ( 6-2*abs(tau-tau_prime)+(tau-tau_prime)**2 - \
            abs(tau - tau_prime)**3 ) ) - math.exp(-t_max) * ( 9*t_max**3 + \
            27*t_max**2 + 54*t_max + 54 - 24*(tau**2+tau_prime**2)*(t_max+1) + \
            12*(tau**2-tau_prime**2)**2*(tau**2+tau_prime**2)*(t_max-1)/t_max**2 + \
            3*(tau**2-tau_prime**2)**4/(8*t_max**4)*(6-2*t_max+t_max**2-t_max**3) ) + \
            ( 30*tau**4 + 20*tau**2*tau_prime**2 + 30*tau_prime**4 - \
            12*(tau**2-tau_prime**2)**2*(tau**2+tau_prime**2) + \
            3*(tau**2-tau_prime**2)**4/8 ) * \
            ( -special.exp1(t_max) ) )

print('Functions B00, B20, f02, f20, f22 are calculated.')    

for tau_p in range(0, height, 1): # I00, I02, I20, I22
    tau = tau_p / (height-1) * optical_radius
    tau += tau_star

    if(tau == tau_star):
        delta_1 = 0
        delta_2 = 0.5 / (height-1) * optical_radius
    elif(tau_p == height-1):
        delta_1 = 0.5 / (height-1) * optical_radius
        delta_2 = 0
    else:
        delta_1 = 0.5 / (height-1) * optical_radius
        delta_2 = 0.5 / (height-1) * optical_radius

    # expression (40)
    I00_array[tau_p] = albedo + albedo/(2*tau)*(special.expn(3, delta_1) - special.expn(3, delta_2)) - \
        albedo/2*(1 - delta_1/tau)*special.expn(2, delta_1) - albedo/2*(1+delta_2/tau)*special.expn(2, delta_2) + \
        albedo/(2*tau)*math.sqrt( (tau+delta_2)**2 - tau_star**2) * \
        special.expn(2, math.sqrt(tau**2 - tau_star**2) + math.sqrt( (tau + delta_2)**2 - tau_star**2) ) + \
        albedo/(2*tau)*special.expn(3, math.sqrt(tau**2 - tau_star**2) + math.sqrt( (tau + delta_2)**2 - tau_star**2) ) - \
        albedo/(2*tau)*math.sqrt( (tau - delta_1)**2 - tau_star**2 ) * \
        special.expn(2, math.sqrt(tau**2 - tau_star**2) + math.sqrt( (tau - delta_1)**2 - tau_star**2) ) - \
        albedo/(2*tau)*special.expn(3, math.sqrt(tau**2 - tau_star**2) + math.sqrt( (tau - delta_1)**2 - tau_star**2) )

    # expression (42)
    I02_array[tau_p] = albedo*( 16/(15*tau) + ( (tau+delta_1)/4 - ( 1/(8*tau) + 3*tau/8 )*delta_1**2 + 3/8*delta_1**3 - \
        3/(32*tau)*delta_1**4 ) * special.expn(2, delta_1) + \
        (1/4 - ( 1/(4*tau) + 3/4*tau )*delta_1 + 9/8*delta_1**2 - 3/8/tau*delta_1**3 ) * special.expn(3, delta_1) + \
        (-1/4/tau - 3/4*tau + 9/4*delta_1 - 9/8/tau*delta_1**2 ) * special.expn(4, delta_1) + \
        9/4*(1-delta_1/tau)*special.expn(5, delta_1) - 9/4/tau*special.expn(6, delta_1) + \
        ( (tau-delta_2)/4 - ( 1/8/tau + 3/8*tau)*delta_2**2 - 3/8*delta_2**3 - 3/32/tau*delta_2**4 ) * special.expn(2, delta_2) - \
        ( 1/4 + (1/4/tau + 3/4*tau)*delta_2 + 9/8*delta_2**2 + 3/8/tau*delta_2**3 ) * special.expn(3, delta_2) - \
        ( 1/4/tau + 3/4*tau + 9/4*delta_2 + 9/8/tau*delta_2**2) * special.expn(4, delta_2) - \
        9/4*(1 + delta_2/tau)*special.expn(5, delta_2) - 9/4/tau*special.expn(6, delta_2) )

    # expression (45)
    I20_array[tau_p] = albedo * ( 3/2/tau**2 + ( 1/8 - delta_1/2/tau + (9/16/tau**2 - 3/16)*delta_1**2 + (3/8/tau - 3/16/tau**3) * \
        delta_1**3 - 15*delta_1**4/64/tau**2 + 3*delta_1**5/64/tau**3) * special.expn(2, delta_1) + \
        ( -1/2/tau + (9/8/tau**2 - 3/8)*delta_1 + (9/8/tau - 9/16/tau**3)*delta_1**2 - \
        15*delta_1**3/16/tau**2 + 15*delta_1**4/64/tau**3) * special.expn(3, delta_1) + \
        ( 9/8/tau**2 - 3/8 + (9/4/tau - 9/8/tau**3)*delta_1 - 45*delta_1**2/16/tau**2 + 15*delta_1**3/16/tau**3 ) * \
        special.expn(4, delta_1) + \
        ( 9/4/tau - 9/8/tau**3 - 45*delta_1/8/tau**2 + 45*delta_1**2/16/tau**3 ) * special.expn(5, delta_1) + \
        ( -45/8/tau**2 + 45*delta_1/8/tau**3 ) * special.expn(6, delta_1) + 45/8/tau**3 * special.expn(7, delta_1) + \
        ( 1/8 + delta_2/2/tau + ( 9/16/tau**2 - 3/16 )*delta_2**2 + ( 3/16/tau**3 - 3/8/tau ) * delta_2**3 - \
        15*delta_2**4/64/tau**2 - 3*delta_2**5/64/tau**3 ) * special.expn(2, delta_2) + \
        ( 1/2/tau + ( 9/8/tau**2 - 3/8 )*delta_2 + ( 9/16/tau**3 - 9/8/tau )*delta_2**2 - \
        15*delta_2**3/16/tau**2 - 15*delta_2**4/64/tau**3 ) * special.expn(3, delta_2) + \
        ( 9/8/tau**2 - 3/8 + ( 9/8/tau**3 - 9/4/tau )*delta_2 - 45*delta_2**2/16/tau**2 - 15*delta_2**3/16/tau**3 ) * \
        special.expn(4, delta_2) + \
        ( 9/8/tau**3 - 9/4/tau - 45*delta_2/8/tau**2 - 45*delta_2**2/16/tau**3) * special.expn(5, delta_2) + \
        ( -45/8/tau**2 - 45*delta_2/8/tau**3 ) * special.expn(6, delta_2) - 45/8/tau**3 * special.expn(7, delta_2) )

    # expression (47)
    I22_array[tau_p] = albedo*( 129/14/tau**3 + 829/210/tau + 7/10*tau + ( -5/8*tau + 5/4*delta_1 + (3/4*tau - 25/16/tau) * \
        delta_1**2 + (15/16/tau**2 - 3/2)*delta_1**3 + ( -15/64/tau**3 + 21/16/tau - 3/64*tau )*delta_1**4 + \
        ( 3/32 - 9/16/tau**2 )*delta_1**5 + ( 3/32/tau**3 - 9/128/tau )*delta_1**6 + \
        3/128/tau**2*delta_1**7 - 3/1024/tau**3*delta_1**8 ) * special.expn(2, delta_1) + \
        ( 5/4 + ( 3/2*tau - 25/8/tau )*delta_1 + ( 45/16/tau**2 - 9/2 )*delta_1**2 + \
        ( -15/16/tau**3 + 21/4/tau - 3/16*tau )*delta_1**3 + \
        ( 15/32 - 45/16/tau**2 )*delta_1**4 + ( 9/16/tau**3 - 27/64/tau )*delta_1**5 + \
        21/128/tau**2*delta_1**6 - 3/128/tau**3*delta_1**7 ) * special.expn(3, delta_1) + \
        ( 3/2*tau - 25/8/tau + ( 45/8/tau**2 - 9 )*delta_1 + \
        ( -45/16/tau**3 + 63/4/tau - 9/16*tau )*delta_1**2 + \
        ( 15/8 - 45/4/tau**2 )*delta_1**3 + ( 45/16/tau**3 - 135/64/tau )*delta_1**4 + \
        63/64/tau**2*delta_1**5 - 21/128/tau**3*delta_1**6 ) * special.expn(4, delta_1) + \
        ( 45/8/tau**2 - 9 + ( -45/8/tau**3 + 63/2/tau - 9/8*tau )*delta_1 + \
        ( 45/8 - 135/4/tau**2 )*delta_1**2 + \
        ( 45/4/tau**3 - 135/16/tau )*delta_1**3 + 315/64/tau**2*delta_1**4 - \
        63/64/tau**3*delta_1**5 ) * special.expn(5, delta_1) + \
        ( -45/8/tau**3 + 63/2/tau - 9/8*tau + ( 45/4 - 135/2/tau**2 )*delta_1 + \
        ( 135/4/tau**3 - 405/16/tau )*delta_1**2 + \
        315/16/tau**2*delta_1**3 - 315/64/tau**3*delta_1**4 ) * special.expn(6, delta_1) + \
        ( 45/4 - 135/2/tau**2 + ( 135/2/tau**3 - 405/8/tau )*delta_1 + 945/16/tau**2*delta_1**2 - \
        315/16/tau**3*delta_1**3 ) * special.expn(7, delta_1) + \
        ( 135/2/tau**3 - 405/8/tau + 945/8/tau**2*delta_1 - 945/16/tau**3*delta_1**2 ) * \
        special.expn(8, delta_1) + \
        945/8*( 1/tau**2 - delta_1/tau**3 ) * special.expn(9, delta_1) - \
        945/8/tau**3 * special.expn(10, delta_1) + \
        ( -5/8*tau - 5/4*delta_2 + ( 3/4*tau - 25/16/tau )*delta_2**2 + ( 3/2 - 15/16/tau**2)*delta_2**3 + \
        ( -15/64/tau**3 + 21/16/tau - 3/64*tau )*delta_2**4 + ( 9/16/tau**2 - 3/32 )*delta_2**5 + \
        ( 3/32/tau**3 - 9/128/tau )*delta_2**6 - \
        3/128/tau**2*delta_2**7 - 3/1024/tau**3*delta_2**8 ) * special.expn(2, delta_2) + \
        ( -5/4 + ( 3/2*tau - 25/8/tau )*delta_2 + ( 9/2 - 45/16/tau**2 )*delta_2**2 + \
        ( -15/16/tau**3 + 21/4/tau - 3/16*tau )*delta_2**3 + \
        ( 45/16/tau**2 - 15/32 )*delta_2**4 + ( 9/16/tau**3 - 27/64/tau )*delta_2**5 - \
        21/128/tau**2*delta_2**6 - 3/128/tau**3*delta_2**7 ) * special.expn(3, delta_2) + \
        ( 3/2*tau - 25/8/tau + ( 9 - 45/8/tau**2)*delta_2 + ( -45/16/tau**3 + \
        63/4/tau - 9/16*tau )*delta_2**2 + \
        ( 45/4/tau**2 - 15/8 )*delta_2**3 + ( 45/16/tau**3 - 135/64/tau )*delta_2**4 - \
        63/64/tau**2*delta_2**5 - 21/128/tau**3*delta_2**6 ) * special.expn(4, delta_2) + \
        ( 9 - 45/8/tau**2 + ( -45/8/tau**3 + 63/2/tau - 9/8*tau )*delta_2 + \
        ( 135/4/tau**2 - 45/8 )*delta_2**2 + \
        ( 45/4/tau**3 - 135/16/tau )*delta_2**3 - 315/64/tau**2*delta_2**4 - \
        63/64/tau**3*delta_2**5) * special.expn(5, delta_2) + \
        ( -45/8/tau**3 + 63/2/tau - 9/8*tau + ( 135/2/tau**2 - 45/4 )*delta_2 + \
        ( 135/4/tau**3 - 405/16/tau )*delta_2**2 - \
        315/16/tau**2*delta_2**3 - 315/64/tau**3*delta_2**4 ) * special.expn(6, delta_2) + \
        ( 135/2/tau**2 - 45/4 + ( 135/2/tau**3 - 405/8/tau )*delta_2 - 945/16/tau**2*delta_2**2 - \
        315/16/tau**3*delta_2**3 ) * special.expn(7, delta_2) + \
        ( 135/2/tau**3 - 405/8/tau - 945/8/tau**2*delta_2 - 945/16/tau**3*delta_2**2 ) *
        special.expn(8, delta_2) - \
        945/8 * (1/tau**2 + delta_2/tau**3) * special.expn(9, delta_2) - \
        945/8/tau**3 * special.expn(10, delta_2) )


print('Functions I00, I02, I20, I22 are calculated.')

# solving system of linear equations to find B0 and B2

twodim_array_linear = [[0 for x in range(width*2)] for y in range(height*2)] # matrix of coefficients for linear equations
B00_B20_array_linear = [0 for y in range(height*2)] # vector of constants for linear equations

for tau_p in range(0, height, 1):
    B00_B20_array_linear[tau_p] = -B00_array[tau_p]
    B00_B20_array_linear[height + tau_p] = -B20_array[tau_p]

B00_B20_array_linear = np.array(B00_B20_array_linear)

for tau_p in range(0, height, 1):
    for tau_prime_p in range(0, width, 1):

        tau = tau_p / (height-1) * optical_radius
        tau += tau_star
        
        if(tau_prime_p == 0 or tau_prime_p == (width-1)):
            int_weight = 1 # first or last element for Simpsons method
        elif tau_prime_p % 2 == 1 :
            int_weight = 4 # odd element for Simpsons method
        else :
            int_weight = 2 # even element for Simpsons method

        # 1/(width-1)/3*optical_radius is a step for Simpsons method

        twodim_array_linear[tau_p][tau_prime_p] = int_weight * K00_array[tau_p][tau_prime_p] * \
                                                   1/(width-1)/3 * optical_radius                # B0 K00
        
        twodim_array_linear[tau_p][width+tau_prime_p] = int_weight * K02_array[tau_p][tau_prime_p] * \
                                                         1/(width-1)/3 * optical_radius          # B0 K02

        twodim_array_linear[height+tau_p][tau_prime_p] = int_weight * K20_array[tau_p][tau_prime_p] * \
                                                          1/(width-1)/3 * optical_radius         # B2 K20

        twodim_array_linear[height+tau_p][width+tau_prime_p] = int_weight * K22_array[tau_p][tau_prime_p] * \
                                                                1/(width-1)/3 * optical_radius   # B2 K22

        if(tau_p == tau_prime_p):
            twodim_array_linear[tau_p][tau_prime_p] = -1 + I00_array[tau_p]                                        # -B0(tau_p)

            twodim_array_linear[tau_p][width+tau_prime_p] = int_weight*1/(width-1)/3*optical_radius*f02_array[tau_p] + \
                                                             I02_array[tau_p]/tau

            twodim_array_linear[height+tau_p][tau_prime_p] = int_weight*1/(width-1)/3*optical_radius*f20_array[tau_p] + \
                                                              I20_array[tau_p]                                                            
            
            twodim_array_linear[height+tau_p][width+tau_prime_p] = -1 + int_weight*1/(width-1)/3*optical_radius*f22_array[tau_p] + \
                                                                    I22_array[tau_p]/tau                            # -B2(tau_p)
            

            
twodim_array_linear = np.array(twodim_array_linear)

B0_B2_array = np.linalg.solve(twodim_array_linear, B00_B20_array_linear)

del twodim_array_linear


# checking B0 and B2 results for physical correctness according to formulas (10),(11)
for tau_p in range(0, height, 1):
    B0_array[tau_p] = B0_B2_array[tau_p]
    B2_array[tau_p] = B0_B2_array[height + tau_p]
    if(B0_array[tau_p] < 0):
        print('Error, B0_array[', tau_p, '] is less than 0')
        print(B0_array[tau_p])
    if(B2_array[tau_p] < -B0_array[tau_p] or B2_array[tau_p] > 0.5*B0_array[tau_p]):
        print('Eror, B2_array[', tau_p, '] = ', B2_array[tau_p] )
        print('B0_array[', tau_p, '] = ', B0_array[tau_p] )

del B0_B2_array

print('Source functions B0 and B2 are calculated.')


for mu_p in range(0, height, 1):
    mu = mu_p / (height-1)

    integral_I_mu_1 = 0
    integral_I_mu_2 = 0
    integral_Q_mu_1 = 0
    integral_Q_mu_2 = 0


    if(mu == 0):
        
        I_array[mu_p] = 0 # integral is equal to zero
        Q_array[mu_p] = 0 # integral is equal to zero

    elif(mu > 0 and mu <= math.sqrt(1 - tau_star**2/full_opt_rad**2) ):

        point_quantity = width # point quantity for integral calculation using Simpsons method

        for point_p in range(0, point_quantity, 1):
            
            if(point_p == 0 or point_p == (point_quantity-1)):
                int_weight = 1 # first or last element for Simpsons method
            elif point_p % 2 == 1 :
                int_weight = 4 # odd element for Simpsons method
            else :
                int_weight = 2 # even element for Simpsons method

            step = ( full_opt_rad - full_opt_rad*math.sqrt(1 - mu**2) ) / (point_quantity - 1)

            tau_prime = ( full_opt_rad - full_opt_rad*math.sqrt(1 - mu**2) ) * (point_p+1)/point_quantity + \
                         full_opt_rad*math.sqrt(1 - mu**2)


            point_for_B0_B2 = (tau_prime-tau_star) * (width-1)/optical_radius
           
            if( int(point_for_B0_B2) < (width-1) ):
                B0_array_value = ( point_for_B0_B2 - int(point_for_B0_B2) ) * ( B0_array[int(point_for_B0_B2)+1] - \
                                 B0_array[int(point_for_B0_B2)] ) + B0_array[int(point_for_B0_B2)]

                B2_array_value = ( point_for_B0_B2 - int(point_for_B0_B2) ) * ( B2_array[int(point_for_B0_B2)+1] - \
                                 B2_array[int(point_for_B0_B2)] ) + B2_array[int(point_for_B0_B2)]
            else:
                B0_array_value = B0_array[width-1]
                B2_array_value = B2_array[width-1]
                
            # formula (62)
            integral_I_mu_1 += int_weight*step/3 * 2 * ( B0_array_value + B2_array_value * ( 1 - 3*full_opt_rad**2 * \
                               (1 - mu**2)/2/tau_prime**2) ) * math.exp(-full_opt_rad*mu) * \
                               (math.exp( math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) + \
                               math.exp( -math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) ) / 2 * \
                               tau_prime/math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) )

            # formula (64)
            integral_Q_mu_1 += int_weight*step/3 * ( -3*full_opt_rad**2 * (1 - mu**2) * B2_array_value * \
                               math.exp(-full_opt_rad*mu) *\
                               (math.exp( math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) + \
                               math.exp( -math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) ) / 2 * \
                               1 / tau_prime / math.sqrt( tau_prime**2 - full_opt_rad**2*(1-mu**2) ) )
                                                     
            
        I_array[mu_p] = integral_I_mu_1
        Q_array[mu_p] = integral_Q_mu_1
        
    elif(mu > math.sqrt(1 - tau_star**2/full_opt_rad**2) ):

        for tau_prime_p in range(0, width, 1):

            tau_prime = tau_prime_p / (width-1) * optical_radius
            tau_prime += tau_star

            step = (full_opt_rad - tau_star)/(width-1)

            if(tau_prime_p == 0 or tau_prime_p == (width-1)):
                int_weight = 1 # first or last element for Simpsons method
            elif tau_prime_p % 2 == 1 :
                int_weight = 4 # odd element for Simpsons method
            else :
                int_weight = 2 # even element for Simpsons method

            # part of formula (63)
            integral_I_mu_2 += int_weight*step/3*( B0_array[tau_prime_p] + B2_array[tau_prime_p] * \
                               ( 1 - 3*full_opt_rad**2*(1-mu**2)/2/tau_prime**2 ) ) * \
                               math.exp(-full_opt_rad*mu + math.sqrt(tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) * \
                               tau_prime / math.sqrt( tau_prime**2 - full_opt_rad**2*(1 - mu**2) )

            # formula (65)
            integral_Q_mu_2 += -3/2*full_opt_rad**2*(1-mu**2) * int_weight*step/3 * \
                               B2_array[tau_prime_p] * \
                               math.exp(-full_opt_rad*mu + math.sqrt(tau_prime**2 - full_opt_rad**2*(1-mu**2) ) ) * \
                               1 / tau_prime / math.sqrt( tau_prime**2 - full_opt_rad**2*(1 - mu**2) )

        # formula (63)
        I_array[mu_p] = integral_I_mu_2 + Luminosity*a**2/(4*math.pi**2*full_opt_rad**2) * \
                        math.exp(-full_opt_rad*mu + math.sqrt(tau_star**2 - full_opt_rad**2*(1-mu**2) ) )

        Q_array[mu_p] = integral_Q_mu_2
        
    else:
        
        print('Error - problems with mu parameter')


print('Stokes parameters I and Q are calculated.')

# to get Stokes I parameter
ro_0 = cloud_radius / dist_star_cent_tel

for mu_p in range(0, height, 1):

    mu = mu_p / (height-1)

    ro = math.sqrt(1 - mu**2) * ro_0

    dist_cloud_center = math.tan(ro) * dist_star_cent_tel # radius from cloud center to create images

    rad_array = dist_cloud_center / cloud_radius *39/40*image_x_dim/2  # radius for 2d array for Stokes vector

    for angle_psi_p in range(0, 3600, 1):
        angle_psi = angle_psi_p / 1800 * math.pi
        data_gnu_I_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = I_array[mu_p] 

        data_gnu_Q_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                                       Q_array[mu_p]*math.cos(2*angle_psi)

        data_gnu_U_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                                       -Q_array[mu_p]*math.sin(2*angle_psi)

# four central pixels are center of the star
data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2 - 1)] = data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_I_2[int(image_y_dim/2 - 1)][int(image_x_dim/2)] = data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_I_2[int(image_y_dim/2 - 1)][int(image_x_dim/2 - 1)] = data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2)]

data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2 - 1)] = data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_Q_2[int(image_y_dim/2 - 1)][int(image_x_dim/2)] = data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_Q_2[int(image_y_dim/2 - 1)][int(image_x_dim/2 - 1)] = data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2)]

data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2 - 1)] = data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_U_2[int(image_y_dim/2 - 1)][int(image_x_dim/2)] = data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2)]
data_gnu_U_2[int(image_y_dim/2 - 1)][int(image_x_dim/2 - 1)] = data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2)]

'''
# to get several Stokes parameters I,Q,U using linear interpolation    
for rad_array in range(0, 10, 1):
    if(data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)] == 0 ): # if one central cloud pixel is zero

        dist_cloud_center = rad_array*cloud_radius /39*40/image_x_dim/2
        ro = math.atan2(dist_cloud_center, dist_star_cent_tel)
        mu = math.sqrt(1 - ro**2/ro_0**2)

        point_mu = mu*(width-1) # index for Stokes vector parameters

        if( int(point_mu) < (width-1) ):
            I_array_val = (point_mu - int(point_mu) ) * ( I_array[int(point_mu) + 1] - \
                          I_array[int(point_mu)] ) + I_array[int(point_mu)]

            Q_array_val = (point_mu - int(point_mu) ) * ( Q_array[int(point_mu) + 1] - \
                            Q_array[int(point_mu)] ) + Q_array[int(point_mu)] 
                            
        else:
            I_array_val = I_array[int(width-1)]
            Q_array_val = Q_array[int(width-1)]

        for angle_psi_p in range(0, 3600, 1):
            angle_psi = angle_psi_p / 1800 * math.pi
            
            data_gnu_I_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = I_array_val

            data_gnu_Q_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                                           Q_array_val*math.cos(2*angle_psi)

            data_gnu_U_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                                           -Q_array_val*math.sin(2*angle_psi)

for rad_array in range(0, 10, 1):
    for angle_psi_p in range(0, 3600, 1):
        angle_psi = angle_psi_p / 1800 * math.pi

        data_gnu_I_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                data_gnu_I_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)]

        data_gnu_Q_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)]*math.cos(2*angle_psi) + \
                data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)]*math.sin(2*angle_psi)

        data_gnu_U_2[int(image_y_dim/2 + rad_array*math.sin(angle_psi))][int(image_x_dim/2 + rad_array*math.cos(angle_psi))] = \
                -data_gnu_Q_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)]*math.sin(2*angle_psi) + \
                data_gnu_U_2[int(image_y_dim/2)][int(image_x_dim/2 + rad_array)]*math.cos(2*angle_psi)
'''

file = open('data_gnu_I_2.txt', 'w')
for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        file.write("{:.4e}".format(data_gnu_I_2[y][x]))
        file.write(' ')
    file.write('\n')
file.close()

file = open('data_gnu_Q_2.txt', 'w')
for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        file.write("{:.4e}".format(data_gnu_Q_2[y][x]))
        file.write(' ')
    file.write('\n')
file.close()

file = open('data_gnu_U_2.txt', 'w')
for y in range(0, image_y_dim, 1):
    for x in range(0, image_x_dim, 1):
        file.write("{:.4e}".format(data_gnu_U_2[y][x]))
        file.write(' ')
    file.write('\n')
file.close()



file = open('file_001_I_array.txt', 'w')
for mu_p in range(0, height, 1):
    file.write("{:.4e}".format(I_array[mu_p]))
    file.write('\n')
file.close()

file = open('file_002_Q_array.txt', 'w')
for mu_p in range(0, height, 1):
    file.write("{:.4e}".format(Q_array[mu_p]))
    file.write('\n')
file.close()



file = open('file_005_B0.txt', 'w')
for y in range(0, height, 1): # tau
    file.write("{:.4e}".format(B0_array[y]))
    file.write('\n')
file.close()

file = open('file_006_B2.txt', 'w')
for y in range(0, height, 1): # tau
    file.write("{:.4e}".format(B2_array[y]))
    file.write('\n')
file.close()


print('Program finished the calculations.')


