/**
 * @file      photon_packet.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Photon_packet class functions.
 */


#include "photon_packet.h"


Photon_packet::Photon_packet()
{
    wave_length = 0; 
    
    stokes_I_param = 0; 
    stokes_q_param = 0; 
    stokes_u_param = 0; 
    stokes_v_param = 0; 
    
    PHI_1 = 0;           
    THET_1 = 0;          

    prev_x = 0;
    prev_y = 0;
    prev_z = 0;

    geometrical_way = 0;

    optical_way_0 = 10;

    optical_way = 0;

    photon_location = 0;

    photon_weight = 0;
    photon_energy = 0;  

    cloud_type = 0;
    concentration_type = 0;

    inner_surface = 0;
    outer_surface = 0;
    inside_cloud = 0;

    stokes_q_tel = 0;
    stokes_u_tel = 0;
    stokes_v_tel = 0;
}



void Photon_packet::set_wave_length(double length)
{
    wave_length = length;
}


double Photon_packet::get_wave_length()
{
    return wave_length;
}


void Photon_packet::set_stokes_I_param(double I_par)
{
    stokes_I_param = I_par;
}


double Photon_packet::get_stokes_I_param()
{
    return stokes_I_param;
}


void Photon_packet::set_stokes_q_param(double q_par)
{
    stokes_q_param = q_par;
}


double Photon_packet::get_stokes_q_param()
{
    return stokes_q_param;
}


void Photon_packet::set_stokes_u_param(double u_par)
{
    stokes_u_param = u_par;
}


double Photon_packet::get_stokes_u_param()
{
    return stokes_u_param;
}


void Photon_packet::set_stokes_v_param(double v_par)
{
    stokes_v_param = v_par;
}


double Photon_packet::get_stokes_v_param()
{
    return stokes_v_param;
}


void Photon_packet::set_optical_way(float opt_leng)
{
    optical_way = opt_leng;
}


float Photon_packet::get_optical_way()
{
    return optical_way;
}


void Photon_packet::set_photon_location(short location)
{
    photon_location = location;
}


short Photon_packet::get_photon_location()
{
    return photon_location;
}


void Photon_packet::set_coord_X_current(float x_cor)
{
    coord_X_current = x_cor;
}


float Photon_packet::get_coord_X_current()
{
    return coord_X_current;
}


void Photon_packet::set_coord_Y_current(float y_cor)
{
    coord_Y_current = y_cor;
}


float Photon_packet::get_coord_Y_current()
{
    return coord_Y_current;
}


void Photon_packet::set_coord_Z_current(float z_cor)
{
    coord_Z_current = z_cor;
}


float Photon_packet::get_coord_Z_current()
{
    return coord_Z_current;
}


void Photon_packet::set_coord_X_previous(float x_coor)
{
    coord_X_previous = x_coor;
}


float Photon_packet::get_coord_X_previous()
{
    return coord_X_previous;
}


void Photon_packet::set_coord_Y_previous(float y_coor)
{
    coord_Y_previous = y_coor;
}


float Photon_packet::get_coord_Y_previous()
{
    return coord_Y_previous;
}


void Photon_packet::set_coord_Z_previous(float z_coor)
{
    coord_Z_previous = z_coor;
}


float Photon_packet::get_coord_Z_previous()
{
    return coord_Z_previous;
}


void Photon_packet::set_photon_weight(double phot_weight)
{
    photon_weight = phot_weight;
}


double Photon_packet::get_photon_weight()
{
    return photon_weight;
}


void Photon_packet::set_photon_energy(double phot_energ)
{
    photon_energy = phot_energ;
}


double Photon_packet::get_photon_energy()
{
    return photon_energy;
}


//void Photon_packet::calc_phot_energ_0()
//{
//    photon_energy = photon_quantity * PLANCK_CONST * c_0 / wave_length;
//}


void Photon_packet::sca_phot_pack(double albedo)
{
    //if( 1 + 0.0001 < sqrt(stokes_q_tel*stokes_q_tel + stokes_u_tel*stokes_u_tel + stokes_v_tel*stokes_v_tel))
    //{
    //    photon_weight = 0;
    //}

    //photon_energy *= albedo;
    photon_weight *= albedo;

    calc_Prob_phot_tel();
    phot_weight_tel = photon_weight * Prob_phot_tel * exp(-1.0*optical_way_peel);
    
    double t_cross = distance_tel - x_tel_vect*coord_X_current - y_tel_vect*coord_Y_current - z_tel_vect*coord_Z_current;

    double x_sph_tel_1 = coord_X_current + x_tel_vect*t_cross;
    double y_sph_tel_1 = coord_Y_current + y_tel_vect*t_cross;
    double z_sph_tel_1 = coord_Z_current + z_tel_vect*t_cross;

    double x_sph_tel_2 = x_sph_tel_1 - ARRAY_X/2;
    double y_sph_tel_2 = y_sph_tel_1 - ARRAY_Y/2;
    double z_sph_tel_2 = z_sph_tel_1;

    double x_sph_tel_3 = x_sph_tel_2*cos(-azimuth_tel_pos) - y_sph_tel_2*sin(-azimuth_tel_pos);
    double y_sph_tel_3 = x_sph_tel_2*sin(-azimuth_tel_pos) + y_sph_tel_2*cos(-azimuth_tel_pos);
    double z_sph_tel_3 = z_sph_tel_2;

    double x_sph_tel_4 = x_sph_tel_3 + ARRAY_X/2;
    double y_sph_tel_4 = y_sph_tel_3 + ARRAY_Y/2;
    double z_sph_tel_4 = z_sph_tel_3;

    double x_sph_tel_5 = x_sph_tel_4 - ARRAY_X/2;
    double y_sph_tel_5 = y_sph_tel_4;
    double z_sph_tel_5 = z_sph_tel_4 - ARRAY_Z/2;

    double x_sph_tel_6 = x_sph_tel_5*cos(-polar_tel_pos) + z_sph_tel_5*sin(-polar_tel_pos);
    double y_sph_tel_6 = y_sph_tel_5;
    double z_sph_tel_6 = -x_sph_tel_5*sin(-polar_tel_pos) + z_sph_tel_5*cos(-polar_tel_pos);

    double x_sph_tel_7 = x_sph_tel_6 + ARRAY_X/2;
    double y_sph_tel_7 = y_sph_tel_6;
    double z_sph_tel_7 = z_sph_tel_6 + ARRAY_Z/2 - distance_tel_0;
    

    //
    matr_point->add_cell_phot_weight((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel);
    
    matr_point->add_cell_q_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_q_tel);

    matr_point->add_cell_u_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_u_tel);

    matr_point->add_cell_v_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_v_tel);
    //

    if(photon_weight < 1.0)
    {
        photon_location = 4;
    }
}


void Photon_packet::first_scatter_phot_pack(double albedo)
{ 
    //if( 1 + 0.0001 < sqrt(stokes_q_tel*stokes_q_tel + stokes_u_tel*stokes_u_tel + stokes_v_tel*stokes_v_tel))
    //{
    //    photon_weight = 0;
    //}

    //photon_energy = photon_energy * albedo * (1.0 - exp(-1.0*optical_way_0));
    photon_weight = photon_weight * albedo * (1.0 - exp(-1.0*optical_way_0));   

    calc_Prob_phot_tel();
    
    phot_weight_tel = photon_weight * Prob_phot_tel * exp(-1.0*optical_way_peel);

    double t_cross = distance_tel - x_tel_vect*coord_X_current - y_tel_vect*coord_Y_current - z_tel_vect*coord_Z_current;

    double x_sph_tel_1 = coord_X_current + x_tel_vect*t_cross;
    double y_sph_tel_1 = coord_Y_current + y_tel_vect*t_cross;
    double z_sph_tel_1 = coord_Z_current + z_tel_vect*t_cross;

    double x_sph_tel_2 = x_sph_tel_1 - ARRAY_X/2;
    double y_sph_tel_2 = y_sph_tel_1 - ARRAY_Y/2;
    double z_sph_tel_2 = z_sph_tel_1;

    double x_sph_tel_3 = x_sph_tel_2*cos(-azimuth_tel_pos) - y_sph_tel_2*sin(-azimuth_tel_pos);
    double y_sph_tel_3 = x_sph_tel_2*sin(-azimuth_tel_pos) + y_sph_tel_2*cos(-azimuth_tel_pos);
    double z_sph_tel_3 = z_sph_tel_2;

    double x_sph_tel_4 = x_sph_tel_3 + ARRAY_X/2;
    double y_sph_tel_4 = y_sph_tel_3 + ARRAY_Y/2;
    double z_sph_tel_4 = z_sph_tel_3;

    double x_sph_tel_5 = x_sph_tel_4 - ARRAY_X/2;
    double y_sph_tel_5 = y_sph_tel_4;
    double z_sph_tel_5 = z_sph_tel_4 - ARRAY_Z/2;

    double x_sph_tel_6 = x_sph_tel_5*cos(-polar_tel_pos) + z_sph_tel_5*sin(-polar_tel_pos);
    double y_sph_tel_6 = y_sph_tel_5;
    double z_sph_tel_6 = -x_sph_tel_5*sin(-polar_tel_pos) + z_sph_tel_5*cos(-polar_tel_pos);

    double x_sph_tel_7 = x_sph_tel_6 + ARRAY_X/2;
    double y_sph_tel_7 = y_sph_tel_6;
    double z_sph_tel_7 = z_sph_tel_6 + ARRAY_Z/2 - distance_tel_0;
    

    //
    matr_point->add_cell_phot_weight((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel);
    
    matr_point->add_cell_q_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_q_tel);

    matr_point->add_cell_u_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_u_tel);

    matr_point->add_cell_v_sto((int) (x_sph_tel_7*CCD_PIXELS), 
                                            (int)(y_sph_tel_7*CCD_PIXELS), phot_weight_tel*stokes_v_tel);
    //
}


void Photon_packet::set_PHI_1(double s_PHI_1)
{
    PHI_1 = s_PHI_1;
} 


double Photon_packet::get_PHI_1()
{
    return PHI_1;
}


void Photon_packet::set_THET_1(double s_THET_1)
{
    THET_1 = s_THET_1; 
}


double Photon_packet::get_THET_1()
{
    return THET_1;
}


void Photon_packet::initialize_packet(emit_phot_param pho_struct)
{
    PHI_1 = pho_struct.s_PHI_1;
    THET_1 = pho_struct.s_THET_1;
    stokes_I_param = pho_struct.stokes_I;
    stokes_q_param = pho_struct.stokes_Q / pho_struct.stokes_I;
    stokes_u_param = pho_struct.stokes_U / pho_struct.stokes_I;
    stokes_v_param = pho_struct.stokes_V / pho_struct.stokes_I;

    photon_weight = pho_struct.phot_weight;
    wave_length = pho_struct.wave_leng;

    coord_X_current = pho_struct.phot_pack_start_X;
    coord_Y_current = pho_struct.phot_pack_start_Y;
    coord_Z_current = pho_struct.phot_pack_start_Z;

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;
}


float Photon_packet::max_optical_way_calc(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z])
{
    //float step_for_start = 0.005;
    //float step_for_start = 0.05;
    float step_for_start = 0.1;
    float t_var_ = step_for_start;
    
    photon_vect vect;

    photon_vect cloud_vect; // vector for cloud, to calculate the radius for dust cloud

    optical_way_0 = 1e-10;

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;

    float start_coord_X = coord_X_current;
    float start_coord_Y = coord_Y_current;
    float start_coord_Z = coord_Z_current;

    float vector_length = 0;

    if ( cloud_type == 1) // cloud - sphere      
    {   while( sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                             + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                             + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2)) < outer_Radius)
        {

            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc(t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;

            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 386 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";
                
                break;
            }

            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL)
            {
                    
                if( arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 3)
                {
                    
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length > inner_Radius)
                    {
                        optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                    }
                }
                else if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 1)
                {
                    // outer suface of cloud
                    //outer_surface++;
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length < outer_Radius)
                    {
                        optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                    }
                }
                else
                {
                    //inside_cloud++;
                    optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                }
            
            }        
        }
    }    
    else if(cloud_type == 2) // cloud assymetrical_1
    {
        bool in_cloud = true;
 
        while( in_cloud )    
        {
            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc(t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;
            
            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 454 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";
                in_cloud = false;
                
                break;
            }

            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL)
            {
                current_surface = arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current]->get_cloud_surface();
               
                switch(current_surface)
                {
                    case 1:
                        // outer surface
                        vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                            + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                            + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                        calc_PHI_THET_assymetrical_1(coord_X_current, coord_Y_current, coord_Z_current);

                        cloud_vect = vector_assymetrical_1();

                        vector_assym_1 = sqrt( (cloud_vect.x_par - ARRAY_X/2)*(cloud_vect.x_par - ARRAY_X/2)
                                                + (cloud_vect.y_par - ARRAY_Y/2)*(cloud_vect.y_par - ARRAY_Y/2) 
                                                + (cloud_vect.z_par - ARRAY_Z/2)*(cloud_vect.z_par - ARRAY_Z/2));

                        if(vector_length < vector_assym_1)
                        {
                            optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                        }
                        else
                        {
                            in_cloud = false;
                        }
                        break;
                    case 2:
                        // between inner and outer surface
                        optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                        break;
                    default: // for inner surfaces
                        // radius compare for different stars is different
                        vector_length = sqrt( (coord_X_current - struct_star_p[current_surface-3]->x_pos)*
                                              (coord_X_current - struct_star_p[current_surface-3]->x_pos)
                                            + (coord_Y_current - struct_star_p[current_surface-3]->y_pos)*
                                              (coord_Y_current - struct_star_p[current_surface-3]->y_pos)
                                            + (coord_Z_current - struct_star_p[current_surface-3]->z_pos)*
                                              (coord_Z_current - struct_star_p[current_surface-3]->z_pos));

                        if(vector_length > struct_star_p[current_surface - 3]->cloud_radius)
                        {
                            optical_way_0 = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],optical_way_0);
                        }
                        break;
                }
            }           
        }
    }

    //cout << "start_coord_X = " << start_coord_X << "\n";
    //cout << "start_coord_Y = " << start_coord_Y << "\n";
    //cout << "start_coord_Z = " << start_coord_Z << "\n";

    coord_X_current = start_coord_X;
    coord_Y_current = start_coord_Y;
    coord_Z_current = start_coord_Z;

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;

    prev_x = 0;
    prev_y = 0;
    prev_z = 0; 

    //cout << " optical_way_0 = " << optical_way_0 << "\n";
    //cout << optical_way_0 - 10.000000 << "\n";
    //cout << "inner_surface int = " << inner_surface << "\n";
    //cout << "inside_cloud int = " << inside_cloud << "\n";
    //cout << "outer_surface int = " << outer_surface << "\n";
    //cout << "float step = " << step_for_start << "\n";
    return optical_way_0;
}


photon_vect Photon_packet::vector_calc( float t_var)
{
    photon_vect vect;

    vect.x_par = t_var*sin(THET_1)*cos(PHI_1) + coord_X_previous;
    vect.y_par = t_var*sin(THET_1)*sin(PHI_1) + coord_Y_previous;
    vect.z_par = t_var*cos(THET_1) + coord_Z_previous;

    return vect;
}


double Photon_packet::calc_geometrical_way(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z])
{    
    // before using this function be sure that THET_1 and PHI_1 are set correctly 
    // after photon_packet emission the class Photon_source generates the THET_1 and PHI_1
    // but further it is essential to use class Angles_montecarlo to generate THET_1 and PHI_1

    double optical_way_integ = 1e-10;

    if(optical_way < 2e-10)
    {
        optical_way = 2e-10;
    }

    geometrical_way = 0;
    
    //float step_for_start = 0.001;
    //float step_for_start = 0.02;
    float step_for_start = 0.1;
    
    if(concentration_type == 3 && sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2)) < inner_Radius + 2 )
    {
        step_for_start = 0.0005;
    }
    else if(concentration_type == 2 && sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                             + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                             + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2)) < inner_Radius + 2)
    {
        step_for_start = 0.0005;
    }

    float t_var_ = step_for_start;
    photon_vect vect; // vector for photon packet

    photon_vect cloud_vect; // vector for cloud, to calculate the radius for dust cloud

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;
    
    float vector_length = 0;

    inner_surface = 0; // how many times geometrical_way was calculated in inner surface
    outer_surface = 0; // how many times geometrical_way was calculated in outer surface

    
    if(cloud_type == 1) // cloud - sphere
    {
        while( (sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2)) < outer_Radius
                                 && optical_way_integ < optical_way) 
                || ( optical_way == 0 && optical_way_integ <= optical_way) ) 
        {
            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc( t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;
 
            /*
            if( step_for_start > 0.0001 && optical_way - optical_way_integ < 0.001)
            {
                step_for_start = 0.0001;
                //cout << "Calculating with a step 0.0001\n";
            }
            */
            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 632 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";

                photon_location = 3; // photon package out of cloud
                
                break;
            }

            // sub_Volume exist
            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL ) 
            {
                if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 3)
                {
                    // inner surface
                    //inner_surface++; 
                 
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length > inner_Radius)
                    {
                        optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                       [(int)coord_Z_current], optical_way_integ);
                    }
               
                }
                else if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 1)
                {
                    // outer surface
                    //outer_surface++;
                
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length < outer_Radius)
                    {
                        optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                           [(int)coord_Z_current], optical_way_integ);
                    }
                    else
                    {
                        photon_location = 3; // photon package out of cloud
                    }
                
                }
                else
                {
                    optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                           [(int)coord_Z_current], optical_way_integ);
                }                   
            }

        }   

    }
    else if(cloud_type == 2) // cloud assymetrical_1
    {
        bool in_cloud = true;

        while( ( in_cloud && optical_way_integ < optical_way) || ( optical_way == 0 && optical_way_integ <= optical_way) )
        {
            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc( t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;
 
            /*
            if( step_for_start > 0.0001 && optical_way - optical_way_integ < 0.001)
            {
                step_for_start = 0.0001;
                //cout << "Calculating with a step 0.0001\n";
            }
            */

            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 722 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";
                
                in_cloud = false;
                photon_location = 3;

                break;
            }

            // sub_Volume exist
            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL ) 
            {
                current_surface = arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface();
                
                switch(current_surface)
                {
                    case 1:
                        // outer surface
                        //outer_surface++;
                        vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                            + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                            + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                        calc_PHI_THET_assymetrical_1(coord_X_current, coord_Y_current, coord_Z_current);

                        cloud_vect = vector_assymetrical_1();

                        vector_assym_1 =  sqrt( (cloud_vect.x_par - ARRAY_X/2)*(cloud_vect.x_par - ARRAY_X/2)
                                                + (cloud_vect.y_par - ARRAY_Y/2)*(cloud_vect.y_par - ARRAY_Y/2) 
                                                + (cloud_vect.z_par - ARRAY_Z/2)*(cloud_vect.z_par - ARRAY_Z/2));

                        if(vector_length < vector_assym_1)
                        {
                            optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                           [(int)coord_Z_current], optical_way_integ);
                        }
                        else
                        {
                            in_cloud = false;
                            photon_location = 3; // photon package out of cloud
                            //cout << "photon package out of cloud \n";
                        }
                        break;
                    case 2:
                        optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                           [(int)coord_Z_current], optical_way_integ);
                        break;
                    default: // for inner surfaces
                        // radius compare for different stars is different
                        vector_length = sqrt( (coord_X_current - struct_star_p[current_surface-3]->x_pos)*
                                              (coord_X_current - struct_star_p[current_surface-3]->x_pos)
                                            + (coord_Y_current - struct_star_p[current_surface-3]->y_pos)*
                                              (coord_Y_current - struct_star_p[current_surface-3]->y_pos)
                                            + (coord_Z_current - struct_star_p[current_surface-3]->z_pos)*
                                              (coord_Z_current - struct_star_p[current_surface-3]->z_pos));

                        if(vector_length > struct_star_p[current_surface - 3]->cloud_radius)
                        {
                            optical_way_integ = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current]
                                                                        [(int)coord_Z_current], optical_way_integ);
                        }
                        break;
                }
                                 
            }
        } 

    }                         
    //cout << "delta x = " << coord_X_current - coord_X_previous << "\n";
    //cout << "delta y = " << coord_Y_current - coord_Y_previous << "\n";
    //cout << "delta z = " << coord_Z_current - coord_Z_previous << "\n";
    
    geometrical_way = sqrt(   ( coord_X_current - coord_X_previous) * ( coord_X_current - coord_X_previous) 
                            + ( coord_Y_current - coord_Y_previous) * ( coord_Y_current - coord_Y_previous) 
                            + ( coord_Z_current - coord_Z_previous) * ( coord_Z_current - coord_Z_previous) );

    //cout << "inner_surface = " << inner_surface << "\n"; 
    //cout << "outer_surface = " << outer_surface << "\n";
    //cout << geometrical_way << "\n";
    return geometrical_way;
}


double Photon_packet::add_optical_integ_trapez(Sub_Volume *pointe, double optical_integ)
{
    if(concentration_type == 1)
    {
        optical_integ += sqrt(  (coord_X_current - prev_x) * (coord_X_current - prev_x)
                              + (coord_Y_current - prev_y) * (coord_Y_current - prev_y)
                              + (coord_Z_current - prev_z) * (coord_Z_current - prev_z))
                            * pointe->get_concentration() * pointe->get_extinction(); 
    }
    else if(concentration_type == 2)
    {
        optical_integ += sqrt(   (coord_X_current - prev_x) * (coord_X_current - prev_x)
                              + (coord_Y_current - prev_y) * (coord_Y_current - prev_y)
                              + (coord_Z_current - prev_z) * (coord_Z_current - prev_z))
                            * pointe->get_concentration() * pointe->get_extinction() 
                            / (( sqrt(   (coord_X_current - ARRAY_X/2) * (coord_X_current - ARRAY_X/2)
                                     +  (coord_Y_current - ARRAY_Y/2) * (coord_Y_current - ARRAY_Y/2)
                                     +  (coord_Z_current - ARRAY_Z/2) * (coord_Z_current - ARRAY_Z/2))
                               + sqrt(   (prev_x - ARRAY_X/2) * (prev_x - ARRAY_X/2)
                                     +  (prev_y - ARRAY_Y/2) * (prev_y - ARRAY_Y/2)
                                     +  (prev_z - ARRAY_Z/2) * (prev_z - ARRAY_Z/2)) ) / 2 );
    }
    else if(concentration_type == 3)
    {
        optical_integ += sqrt(  (coord_X_current - prev_x) * (coord_X_current - prev_x)
                              + (coord_Y_current - prev_y) * (coord_Y_current - prev_y)
                              + (coord_Z_current - prev_z) * (coord_Z_current - prev_z))
                            * pointe->get_concentration() * pointe->get_extinction() 
                            / (( sqrt(   (coord_X_current - ARRAY_X/2) * (coord_X_current - ARRAY_X/2)
                                     +  (coord_Y_current - ARRAY_Y/2) * (coord_Y_current - ARRAY_Y/2)
                                     +  (coord_Z_current - ARRAY_Z/2) * (coord_Z_current - ARRAY_Z/2))
                               + sqrt(   (prev_x - ARRAY_X/2) * (prev_x - ARRAY_X/2)
                                     +  (prev_y - ARRAY_Y/2) * (prev_y - ARRAY_Y/2)
                                     +  (prev_z - ARRAY_Z/2) * (prev_z - ARRAY_Z/2)) ) / 2 )
                            / (( sqrt(   (coord_X_current - ARRAY_X/2) * (coord_X_current - ARRAY_X/2)
                                     +  (coord_Y_current - ARRAY_Y/2) * (coord_Y_current - ARRAY_Y/2)
                                     +  (coord_Z_current - ARRAY_Z/2) * (coord_Z_current - ARRAY_Z/2))
                               + sqrt(   (prev_x - ARRAY_X/2) * (prev_x - ARRAY_X/2)
                                     +  (prev_y - ARRAY_Y/2) * (prev_y - ARRAY_Y/2)
                                     +  (prev_z - ARRAY_Z/2) * (prev_z - ARRAY_Z/2)) ) / 2 );
    }
    else
    {
        cout << "you haven't selected concetration_type!\n";
        return -100;
    }

    return optical_integ;
}



float Photon_packet::gener_first_optical_way()
{
    double xi = ((rand() % 1000001) / 1000000.0);
    //double xi = ((rand() % 1001) / 1000.0);

    if( xi < 1e-100)
    {
        xi = 1e-100;
    }
    //cout << "xi = " << xi << "\n";
    //cout << "optical_way_0 = " << optical_way_0 << "\n";
    
    optical_way = -1.0*log(1 - xi*(1.0-exp(-1.0*optical_way_0)));
    //optical_way = -1.0 * log10(xi);

    //cout << "876 optical_way = " << optical_way << "\n";

    return optical_way;
}



float Photon_packet::gener_optical_way()
{
    double xi = ((rand() % 1000001) / 1000000.0);
    //double xi = ((rand() % 1001) / 1000.0);

    if( xi < 1e-100)
    {
        xi = 1e-100;
    }

    //cout << "xi = " << xi << "\n";
    optical_way = -1.0 * log(xi);
    //cout << "optical_way = " << optical_way << "\n";
    
    while(xi < 0.001)
    {
        xi = ((rand() % 1000001) / 1000000.0);

        //cout << "xi = " << xi << "\n";
        if( xi < 1e-100)
        {
            xi = 1e-100;
        }
        
        optical_way += -1.0 * log(xi);
        //cout << "optical_way = " << optical_way << "\n";
    }
    
    //cout << "911 optical_way = " << optical_way << "\n";
    
    return optical_way;
}



void Photon_packet::set_cloud_type(short cloud_typ)
{
    cloud_type = cloud_typ;
}


void Photon_packet::set_concentration_type(short concent_type)
{
    concentration_type = concent_type;
}



void Photon_packet::check_opt_way_simps(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z])
{
    double optical_way_check = 0;
    int iter_number = 0;
    double step = 0.001;
    double t_var_ = step;

    double step_integral = 0;

    double s_1 = 0;
    double s_2 = 0;

    photon_vect vector;
    
    double r_min = 0;

    double r_current = 0;

    double s_current = 0;

    double r_1 = sqrt((coord_X_previous - ARRAY_X/2)*(coord_X_previous - ARRAY_X/2)
                           +(coord_Y_previous - ARRAY_Y/2)*(coord_Y_previous - ARRAY_Y/2)
                           +(coord_Z_previous - ARRAY_Z/2)*(coord_Z_previous - ARRAY_Z/2));

    double r_2 = sqrt((coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                       + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2)
                       + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

    calc_PHI_THET_central(coord_X_previous, coord_Y_previous, coord_Z_previous);
    double bigThet = acos( cos(THET_0)*cos(THET_1) + sin(THET_0)*sin(THET_1)*cos(PHI_0 - PHI_1) );
    
    if(concentration_type != 1 && concentration_type != 2 && concentration_type != 3 )
    {
        cout << "Please set the correct concentration type. \n";
    }

    if(bigThet > 0 && bigThet <= PI/2)
    {
        while(r_current <= r_2)
        {   
            iter_number++;        
        
            vector = vector_calc(t_var_);
            r_current = sqrt((vector.x_par - ARRAY_X/2)*(vector.x_par - ARRAY_X/2)
                               +(vector.y_par - ARRAY_Y/2)*(vector.y_par - ARRAY_Y/2)
                               +(vector.z_par - ARRAY_Z/2)*(vector.z_par - ARRAY_Z/2));
            t_var_ += step;
        }

        if(iter_number % 2 == 1)
        {
            iter_number++; 
        }
  
        step_integral = (r_2 - r_1) / iter_number;

        r_current = r_1;
    
        for(int i = 0; i <= iter_number; i++)
        {
            if(i == 0 || i == iter_number ) // first or last
            {
                if(concentration_type == 1)
                {
                    optical_way_check += r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 2)
                {
                    optical_way_check += 1 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 3)
                {
                    optical_way_check += 1 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
            }
            else if(i % 2 == 0) // even
            {
                if(concentration_type == 1)
                {
                    optical_way_check += 2 * r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 2)
                {
                    optical_way_check += 2 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 3)
                {
                    optical_way_check += 2 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
            }
            else // odd
            {
                if(concentration_type == 1)
                {
                    optical_way_check += 4 * r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 2)
                {
                    optical_way_check += 4 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
                else if(concentration_type == 3)
                {
                    optical_way_check += 4 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                }
            }

            r_current += step_integral;
        }

        optical_way_check *= arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_concentration()
                                 * arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_extinction()
                                 * step_integral/3;
    }
    else if(bigThet > PI/2)
    {
        r_min = r_2;

        vector = vector_calc(t_var_);

        while(!(fabs(vector.x_par - coord_X_current) < 0.002 && fabs(vector.y_par - coord_Y_current) < 0.002 
                && fabs(vector.z_par - coord_Z_current) < 0.002)  )
        {   
            iter_number++;        
        
            vector = vector_calc(t_var_);
            r_current = sqrt((vector.x_par - ARRAY_X/2)*(vector.x_par - ARRAY_X/2)
                               +(vector.y_par - ARRAY_Y/2)*(vector.y_par - ARRAY_Y/2)
                               +(vector.z_par - ARRAY_Z/2)*(vector.z_par - ARRAY_Z/2));
            t_var_ += step;

            if(r_min > r_current)
            {
                r_min = r_current;
            }
        }

        if(iter_number % 2 == 1)
        {
            iter_number++; 
        }
        
        if(r_min < r_2) // need to calculate two integrals  bigThet > PI/2 && r_2 > r_min
        {
            //cout << "bigThet > PI/2 && r_2 > r_min" << "\n";
            
            if(r_min < inner_Radius)
            {
                cout << "Photon packet crossing section without particles \n"
                     << "r_min < inner_radius       r_min = " << r_min << "\n";
            }
            else // r_min > inner_Radius
            {
                s_1 = -sqrt(r_1*r_1 - r_min*r_min);
                s_2 = sqrt(r_2*r_2 - r_min*r_min);

                step_integral = (s_2 - s_1) / iter_number;
                s_current = s_1;

                for(int i = 0; i <= iter_number; i++)
                {
                    if(i==0 || i == iter_number ) // first or last
                    {
                        if(concentration_type == 1)
                        {
                            optical_way_check += 1;
                        }
                        else if(concentration_type == 2)
                        {
                            optical_way_check += 1 / sqrt(r_min*r_min + s_current*s_current);
                        }
                        else if(concentration_type == 3)
                        {
                            optical_way_check += 1 / (r_min*r_min + s_current*s_current);
                        }
                    }
                    else if(i % 2 == 0) // even
                    {
                        if(concentration_type == 1)
                        {
                            optical_way_check += 2;
                        }
                        else if(concentration_type == 2)
                        {
                            optical_way_check += 2 / sqrt(r_min*r_min + s_current*s_current);
                        }
                        else if(concentration_type == 3)
                        {
                            optical_way_check += 2 / (r_min*r_min + s_current*s_current);
                        }
                    }
                    else // odd
                    {
                        if(concentration_type == 1)
                        {
                            optical_way_check += 4;
                        }
                        else if(concentration_type == 2)
                        {
                            optical_way_check += 4 / sqrt(r_min*r_min + s_current*s_current);
                        }
                        else if(concentration_type == 3)
                        {
                            optical_way_check += 4 / (r_min*r_min + s_current*s_current);
                        }
                    }

                    s_current += step_integral;
                }

                optical_way_check *= arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_concentration()
                                     * arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_extinction()
                                     * step_integral / 3;
                
            }
        }
        else // need to calculate only one integral  bigThet > PI/2 && r_2 == r_min  (r_2 is the least radius)
        {

            //cout << "bigThet > PI/2 && r_2 == r_min" << "\n";
            //cout << "Number of iteration = " << iter_number << "\n";
  
            step_integral = (r_1 - r_2) / iter_number;

            r_current = r_2;
    
            for(int i = 0; i <= iter_number; i++)
            {
                if(i==0 || i == iter_number ) // first or last
                {
                    if(concentration_type == 1)
                    {
                        optical_way_check += r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 2)
                    {
                        optical_way_check += 1 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 3)
                    {
                        optical_way_check += 1 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                }
                else if(i % 2 == 0) // even
                {
                    if(concentration_type == 1)
                    {
                        optical_way_check += 2*r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 2)
                    {
                        optical_way_check += 2 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 3)
                    {
                        optical_way_check += 2 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                }
                else // odd
                {
                    if(concentration_type == 1)
                    {
                        optical_way_check += 4*r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 2)
                    {
                        optical_way_check += 4 / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                    else if(concentration_type == 3)
                    {
                        optical_way_check += 4 / r_current / sqrt(r_current*r_current - r_1*r_1*sin(bigThet)*sin(bigThet));
                    }
                }

                r_current += step_integral;
            }

            optical_way_check *= arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_concentration()
                                 * arr[(int)vector.x_par][(int)vector.y_par][(int)vector.z_par]->get_extinction()
                                 * step_integral/3;
        }    
    }
    
    cout << "simps Recalculated optical way = " << fixed << setprecision(10) << optical_way_check << "\n";
    //cout << "simps x_coordinate = " << vector.x_par << "\n";
    //cout << "simps y_coordinate = " << vector.y_par << "\n";
    //cout << "simps z_coordinate = " << vector.z_par << "\n";
    //cout << "simps Number of iteration = " << iter_number << "\n";
    //cout << "r_1 = " << r_1 << "      r_2 = " << r_2 << "      r_min = " << r_min <<  "\n";  
}



void Photon_packet::check_opt_way_analytical(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z])
{
    double r_1 = sqrt((coord_X_previous - ARRAY_X/2)*(coord_X_previous - ARRAY_X/2)
                       + (coord_Y_previous - ARRAY_Y/2)*(coord_Y_previous - ARRAY_Y/2)
                       + (coord_Z_previous - ARRAY_Z/2)*(coord_Z_previous - ARRAY_Z/2));

    double r_2 = sqrt((coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                       + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2)
                       + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));  

    calc_PHI_THET_central(coord_X_previous, coord_Y_previous, coord_Z_previous);   

    double optical_way_check = 0;

    double bigThet = acos(cos(THET_0)*cos(THET_1) + sin(THET_0)*sin(THET_1)*cos(PHI_0 - PHI_1));
    cout << "bigThet = " << bigThet << "\n";
    
    double step = 0.001;
    double t_var_ = step;
    
    double r_current = 0;
    double r_min = r_2;
    
    photon_vect vector = vector_calc(t_var_);
    
    while(!(fabs(vector.x_par - coord_X_current) < 0.002 && fabs(vector.y_par - coord_Y_current) < 0.002 
            && fabs(vector.z_par - coord_Z_current) < 0.002)  )
    {          
        vector = vector_calc(t_var_);
        r_current = sqrt((vector.x_par - ARRAY_X/2)*(vector.x_par - ARRAY_X/2)
                            +(vector.y_par - ARRAY_Y/2)*(vector.y_par - ARRAY_Y/2)
                            +(vector.z_par - ARRAY_Z/2)*(vector.z_par - ARRAY_Z/2));
        t_var_ += step;

        if(r_min > r_current)
        {
            r_min = r_current;
        }
    }

    cout << "r_1 = " << r_1 << "        r_2 = " << r_2 << "      r_min = " << r_min << "\n";

    if(r_min < inner_Radius)
    {
        cout << "Error r_min < inner_Radius, photon_packet crossing section without particles. \n";
    }
    
    if(concentration_type == 1)
    {
        cout << "Concentration_type = 1,   exinction*concentration = const \n";       

        if( bigThet > PI/2 && (r_2 - r_min) > 0.01 ) // bigThet > PI/2 && r_2 != r_min  two integrals
        {
            //cout << "(r_2 - r_min) > 0.01 \n"; 
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                     * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                     * (sqrt(r_1*r_1 - r_min*r_min) + sqrt(r_2*r_2 - r_min*r_min));  
        }
        else if(r_2 > r_1) // bigThet > 0 && bigThet < PI/2    one integral
        {
            //cout << "r_2 > r_1 \n";
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                     * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                     * (sqrt(r_2*r_2 - r_1*r_1*sin(bigThet)*sin(bigThet)) - r_1*cos(bigThet));
        }
        else if(r_2 < r_1) // one integral
        {
            //cout << "r_2 < r_1 \n";
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                     * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                     * (sqrt(r_1*r_1 - r_1*r_1*sin(bigThet)*sin(bigThet)) 
                                        - sqrt(r_2*r_2 - r_1*r_1*sin(bigThet)*sin(bigThet)) );
        }                
    }
    else if(concentration_type == 2)
    {
        cout << "Concentration_type = 2,   exinction*concentration / radius \n";
        
        if(bigThet > PI/2 && (r_2 - r_min) > 0.01) // calculating two integrals
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * log( (1 / sin(bigThet) + sqrt( 1 / (sin(bigThet)*sin(bigThet)) - 1))
                                        / (r_min / (r_1*sin(bigThet)) ));


            optical_way_check += arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * log( (r_2 / (r_1*sin(bigThet)) + sqrt(r_2*r_2 / (r_1*r_1*sin(bigThet)*sin(bigThet)) - 1) )
                                         / (r_min / (r_1*sin(bigThet))  ));      
        }
        else if(r_1 > r_2) // calculating one integral
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * log( (1 / sin(bigThet) + sqrt( 1 / (sin(bigThet)*sin(bigThet)) - 1))
                                        / (r_2 / (r_1*sin(bigThet)) + sqrt(r_2*r_2 / (r_1*r_1*sin(bigThet)*sin(bigThet)) - 1) )); 
        }
        else if(r_1 < r_2) // calculating one integral
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * log( (r_2 / (r_1*sin(bigThet)) + sqrt(r_2*r_2 / (r_1*r_1*sin(bigThet)*sin(bigThet)) - 1) )
                                         / (1 / sin(bigThet) + sqrt( 1 / (sin(bigThet)*sin(bigThet)) - 1 ) ) );
        }                            
    }  
    else if(concentration_type == 3)
    {
        cout << "Concentration_type = 3,   exinction*concentration / radius^2 \n";

        if(bigThet > PI/2 && (r_2 - r_min) > 0.01) // calculating two integrals
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * 2 / (r_1*sin(bigThet)) * atan(sqrt((1/sin(bigThet)-1) / (1/sin(bigThet)+1) ));
            
            optical_way_check += arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * 2 / (r_1*sin(bigThet)) * atan(sqrt( (r_2 / (r_1*sin(bigThet)) - 1) / (r_2 / (r_1*sin(bigThet)) + 1)));
        }
        else if(r_1 > r_2) // calculating one integral
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * 2 / (r_1*sin(bigThet)) * ( atan(sqrt((1/sin(bigThet)-1) / (1/sin(bigThet)+1) )) 
                                     - atan(sqrt( (r_2 / (r_1*sin(bigThet)) - 1) / (r_2 / (r_1*sin(bigThet)) + 1) ))) ; 
        }
        else if(r_1 < r_2) // calculating one integral
        {
            optical_way_check = arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_concentration()
                                    * arr[(int)coord_X_previous][(int)coord_Y_previous][(int)coord_Z_previous]->get_extinction()
                                    * 2 / (r_1*sin(bigThet)) * ( atan(sqrt( (r_2 / (r_1*sin(bigThet)) - 1) / (r_2 / (r_1*sin(bigThet)) + 1)))
                                      - atan(sqrt( (1/sin(bigThet) - 1) / (1/sin(bigThet) + 1) )));
        }                            
    }  
    else if(concentration_type == 0)
    {
        cout << "Error you didn't selected the concentration type. \n";
    }    

    cout << "analytical Recalculated optical way = " << fixed << setprecision(10)
         << optical_way_check << "\n";
    
}



void Photon_packet::calc_PHI_THET_central(double x_1, double y_1, double z_1)
{
    double radius = sqrt((x_1 - ARRAY_X/2)*(x_1 - ARRAY_X/2) + (y_1 - ARRAY_Y/2)*(y_1 - ARRAY_Y/2)
                          +(z_1 - ARRAY_Z/2)*(z_1 - ARRAY_Z/2));
    
    THET_0 = acos((z_1-ARRAY_Z/2) / radius);
    PHI_0 = atan2(y_1 - ARRAY_Y/2,  x_1 - ARRAY_X/2);

    if (PHI_0 < 0)
    {
        PHI_0 += 2*PI;
    }    
}


void Photon_packet::set_CCD_matr_point(CCD_matr *matr_p)
{
    matr_point = matr_p;
}


photon_vect Photon_packet::vector_assymetrical_1()
{
    photon_vect vect;

    vect.x_par = outer_Radius*cos(THET_0_assym_1) + ARRAY_X/2;
    vect.y_par = outer_Radius*(sin(PHI_0_assym_1)*sin(THET_0_assym_1)-sin(THET_0_assym_1/3)) + ARRAY_Y*3/4;
    vect.z_par = outer_Radius*(cos(PHI_0_assym_1)*sin(THET_0_assym_1)+cos(THET_0_assym_1/2)) + ARRAY_Z/4;

    return vect;
}


void Photon_packet::set_outer_inner_radius(short outer_Rad, short inner_Rad)
{
    outer_Radius = outer_Rad; // the outer radius of the cloud
    inner_Radius = inner_Rad; // the inner radius of the cloud
}


void Photon_packet::calc_PHI_THET_assymetrical_1(double x_1, double y_1, double z_1)
{
    calcul_rad = 20; // calculated radius
    
    if(fabs((x_1 - ARRAY_X/2) / calcul_rad) > 1.0)
    {
        calcul_rad = fabs(x_1 - ARRAY_X/2); // cos is in interval [-1;1]
    }

    THET_0_assym_1 = acos((x_1 - ARRAY_X/2) / calcul_rad);              

    PHI_0_assym_1 = atan2( y_1 - ARRAY_Y*3/4 + calcul_rad*sin(THET_0_assym_1 / 3) , 
                                        z_1 - ARRAY_Z/4 - calcul_rad*cos(THET_0_assym_1 / 2) );                  

    if(PHI_0_assym_1 < 0)
    {
        PHI_0_assym_1 += 2*PI;
    }

    x_par_assym = calcul_rad*cos(THET_0_assym_1) + ARRAY_X/2;
    y_par_assym = calcul_rad*(sin(PHI_0_assym_1)*sin(THET_0_assym_1)-sin(THET_0_assym_1 / 3)) + ARRAY_Y*3/4;
    z_par_assym = calcul_rad*(cos(PHI_0_assym_1)*sin(THET_0_assym_1)+cos(THET_0_assym_1 / 2)) + ARRAY_Z/4;

    calc_assym_1_rad_angles(0.01, x_1, y_1, z_1);
    calc_assym_1_rad_angles(0.001, x_1, y_1, z_1);

}


void Photon_packet::calc_assym_1_rad_angles(double precision, double x_1, double y_1, double z_1)
{
    if(fabs(x_1 - x_par_assym) > 200*precision || fabs(y_1 - y_par_assym) > 200*precision || fabs(z_1 - z_par_assym) > 200*precision )
    {
        calcul_rad -= precision*10;
    }

    while(fabs(x_1 - x_par_assym) > 200*precision || fabs(y_1 - y_par_assym) > 200*precision || fabs(z_1 - z_par_assym) > 200*precision
          || PHI_0_assym_1 != PHI_0_assym_1) // check if the number is not a NaN
    {
        calcul_rad += precision;

        THET_0_assym_1 = acos((x_1 - ARRAY_X/2) / calcul_rad);
 
        PHI_0_assym_1 = atan2( y_1 - ARRAY_Y*3/4 + calcul_rad*sin(THET_0_assym_1 / 3) , 
                                z_1 - ARRAY_Z/4 - calcul_rad*cos(THET_0_assym_1 / 2) );

        if(PHI_0_assym_1 < 0)
        {
            PHI_0_assym_1 += 2*PI;
        }

        x_par_assym = calcul_rad*cos(THET_0_assym_1) + ARRAY_X/2;
        y_par_assym = calcul_rad*(sin(PHI_0_assym_1)*sin(THET_0_assym_1)-sin(THET_0_assym_1 / 3)) + ARRAY_Y*3/4;
        z_par_assym = calcul_rad*(cos(PHI_0_assym_1)*sin(THET_0_assym_1)+cos(THET_0_assym_1 / 2)) + ARRAY_Z/4;

        if(calcul_rad > 34)
        {
            cout << "error in assym_1 radius calculation\n";
            cout << "x_1 = " << x_1 << "    y_1 = " << y_1 << "    z_1 = " << z_1 << "\n";
            break;
        }
    }
}


void Photon_packet::initialize_struct_star_p(star_coordinates *star_p[star_quantity])
{
   struct_star_p = star_p;
/*
   int i = 0;
   while(struct_star_p[i] != NULL)
   {
       cout << "struct_star_p[" << i << "]->cloud_radius = " << struct_star_p[i]->cloud_radius << "\n";
       i++; 
   }
*/
}



void Photon_packet::set_star_numb(short star_n)
{
    star_numb = star_n;
}



void Photon_packet::set_tel_parameters(tel_param tel_p)
{
    tel_solid_angle = tel_p.tel_solid_angle;  
    
    tel_X_coord = tel_p.tel_X_coor;       
    tel_Y_coord = tel_p.tel_Y_coor;       
    tel_Z_coord = tel_p.tel_Z_coor;   

    polar_tel_pos = tel_p.teles_polar_angle;     // Telescope position polar angle
    azimuth_tel_pos = tel_p.teles_azimuth_angle; // Telescope position azimuthal angle

    x_tel_vect = tel_p.tel_x_vect;      // x vector of telescope's XY plane position
    y_tel_vect = tel_p.tel_y_vect;      // y vector of telescope's XY plane position
    z_tel_vect = tel_p.tel_z_vect;      // z vector of telescope's XY plane position

    distance_tel_0 = tel_p.distan_cloud_sys_tel_sys;  // distance between cloud system and telescope system

    // distance after telescopes XY plane rotations (rotations was made by telescope class) 
    distance_tel = tel_p.tel_dist_for_calculations;       
}



void Photon_packet::calc_Thet_peel()
{
    double val_for_acos = cos(THET_1)*cos(thet_tel)+sin(THET_1)*sin(thet_tel)*cos(phi_tel - PHI_1);
    
    if(val_for_acos > 1)
    {
        val_for_acos = 1;
    }
    else if(val_for_acos < -1)
    {
        val_for_acos = -1;
    }

    Thet_peel = acos(val_for_acos);

    if(Thet_peel < 1e-10)
    {
        Thet_peel = 1e-10;
    }
    else if(PI - Thet_peel < 1e-10)
    {
        Thet_peel = PI - 1e-10;
    }

    //cout << "Thet_peel = " << Thet_peel << " = " << Thet_peel/PI*180 << "\n";
    
}



void Photon_packet::set_F11_value(float F11_val)
{
    F11_value = F11_val;
    
    //cout << "F11_value = " << F11_value << "\n";
}


void Photon_packet::set_F12_value(float F12_val)
{
    F12_value = F12_val;
}


void Photon_packet::calc_Prob_phot_tel()
{
    //Prob_phot_tel = F11_value*tel_solid_angle / (4*PI); // for unpolarized light

    Prob_phot_tel = (F11_value + F12_value*cos_2psi_tel*stokes_q_param - F12_value*sin_2psi_tel*stokes_u_param)
                    *tel_solid_angle / (4*PI);         // for polarized light
                                                       // stokes_q_param and stokes_u_param are parameters
                                                       // of incident photon packet before scattering
}



float Photon_packet::optical_way_peel_calc(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z])
{
    //float step_for_start = 0.005;
    //float step_for_start = 0.05;
    float step_for_start = 0.1;
    float t_var_ = step_for_start;
    
    //calc_spher_coor_cur_tel(); // calculates spherical coordinates rad_tel,thet_tel,phi_tel

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;

    float start_coord_X = coord_X_current;
    float start_coord_Y = coord_Y_current;
    float start_coord_Z = coord_Z_current;

    photon_vect vect;
    photon_vect cloud_vect; // vector for cloud, to calculate the radius for dust cloud

    optical_way_peel = 1e-30;

    float vector_length = 0;

    if ( cloud_type == 1) // cloud - sphere      
    {   while( sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                             + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                             + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2)) < outer_Radius)
        {

            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc_peel(t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;


            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 1609 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";

                break;
            }


            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL)
            {
                    
                if( arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 3)
                {
                    
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length > inner_Radius)
                    {
                        optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                                    optical_way_peel);
                    }
                }
                else if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] -> get_cloud_surface() == 1)
                {
                    // outer suface of cloud
                    //outer_surface++;
                    vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                 + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                 + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                    if(vector_length < outer_Radius)
                    {
                        optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                                    optical_way_peel);
                    }
                }
                else
                {
                    //inside_cloud++;
                    optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                                optical_way_peel);
                }
            
            }        
        }
    }    
    else if(cloud_type == 2) // cloud assymetrical_1
    {
        bool in_cloud = true;
 
        while( in_cloud )    
        {
            prev_x = coord_X_current;
            prev_y = coord_Y_current;
            prev_z = coord_Z_current;

            vect = vector_calc_peel(t_var_);

            coord_X_current = vect.x_par;
            coord_Y_current = vect.y_par;
            coord_Z_current = vect.z_par;

            t_var_ += step_for_start;


            if((int)coord_X_current < 0 || (int)coord_X_current > ARRAY_X-1 || (int)coord_Y_current < 0
                || (int)coord_Y_current > ARRAY_Y-1  || (int)coord_Z_current < 0  || (int)coord_Y_current > ARRAY_Z-1 )
            {
                cout << "Out of array size line 1682 \n";
                cout << "(int)x = " << (int)coord_X_current << "   (int)y = " << (int)coord_Y_current
                     << "   (int)z = " << (int)coord_Z_current << "\n";

                cout << "x = " << coord_X_current << "   y = " << coord_Y_current
                     << "   z = " << coord_Z_current << "\n";
                
                in_cloud = false;

                break;
            }


            if(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current] != NULL)
            {
                current_surface = arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current]->get_cloud_surface();
               
               switch(current_surface)
               {
                    case 1:
                        // outer surface
                        vector_length = sqrt( (coord_X_current - ARRAY_X/2)*(coord_X_current - ARRAY_X/2)
                                            + (coord_Y_current - ARRAY_Y/2)*(coord_Y_current - ARRAY_Y/2) 
                                            + (coord_Z_current - ARRAY_Z/2)*(coord_Z_current - ARRAY_Z/2));

                        calc_PHI_THET_assymetrical_1(coord_X_current, coord_Y_current, coord_Z_current);

                        cloud_vect = vector_assymetrical_1();

                        vector_assym_1 = sqrt( (cloud_vect.x_par - ARRAY_X/2)*(cloud_vect.x_par - ARRAY_X/2)
                                                + (cloud_vect.y_par - ARRAY_Y/2)*(cloud_vect.y_par - ARRAY_Y/2) 
                                                + (cloud_vect.z_par - ARRAY_Z/2)*(cloud_vect.z_par - ARRAY_Z/2));

                        if(vector_length < vector_assym_1)
                        {
                            optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                                        optical_way_peel);
                        }
                        else
                        {
                            in_cloud = false;
                        }
                        break;
                    case 2:
                        // between inner and outer surface
                        optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                                    optical_way_peel);
                        break;
                    default: // for inner surfaces
                        // radius compare for different stars is different
                        vector_length = sqrt( (coord_X_current - struct_star_p[current_surface-3]->x_pos)*
                                              (coord_X_current - struct_star_p[current_surface-3]->x_pos)
                                            + (coord_Y_current - struct_star_p[current_surface-3]->y_pos)*
                                              (coord_Y_current - struct_star_p[current_surface-3]->y_pos)
                                            + (coord_Z_current - struct_star_p[current_surface-3]->z_pos)*
                                              (coord_Z_current - struct_star_p[current_surface-3]->z_pos));

                        if(vector_length > struct_star_p[current_surface - 3]->cloud_radius)
                        {
                            optical_way_peel = add_optical_integ_trapez(arr[(int)coord_X_current][(int)coord_Y_current][(int)coord_Z_current],
                                                optical_way_peel);
                        }
                        break;
               }
            }           
        }
    }

    
    //cout << "start_coord_X = " << start_coord_X << "\n";
    //cout << "start_coord_Y = " << start_coord_Y << "\n";
    //cout << "start_coord_Z = " << start_coord_Z << "\n";

    //cout << "coord_X_current = " << coord_X_current << "\n";
    //cout << "coord_Y_current = " << coord_Y_current << "\n";
    //cout << "coord_Z_current = " << coord_Z_current << "\n";

    coord_X_current = start_coord_X;
    coord_Y_current = start_coord_Y;
    coord_Z_current = start_coord_Z;

    coord_X_previous = coord_X_current;
    coord_Y_previous = coord_Y_current;
    coord_Z_previous = coord_Z_current;

    prev_x = 0;
    prev_y = 0;
    prev_z = 0; 

    //cout << " optical_way_peel = " << optical_way_peel << "\n";
    //cout << "inner_surface int = " << inner_surface << "\n";
    //cout << "inside_cloud int = " << inside_cloud << "\n";
    //cout << "outer_surface int = " << outer_surface << "\n";
    //cout << "float step = " << step_for_start << "\n";
    return optical_way_peel;
}



void Photon_packet::calc_spher_coor_cur_tel()
{
    rad_tel = sqrt( ((float)tel_X_coord-coord_X_current)*((float)tel_X_coord-coord_X_current)
                                +((float)tel_Y_coord-coord_Y_current)*((float)tel_Y_coord-coord_Y_current)
                                +((float)tel_Z_coord-coord_Z_current)*((float)tel_Z_coord-coord_Z_current) );

    thet_tel = acos( ((float)tel_Z_coord - coord_Z_current) / rad_tel);   
    phi_tel = atan2( (float)tel_Y_coord - coord_Y_current, (float)tel_X_coord - coord_X_current);

    if(phi_tel < 0)
    {
        phi_tel += 2*PI;
    }
}



photon_vect Photon_packet::vector_calc_peel( float t_var)
{
    photon_vect vect;

    vect.x_par = t_var*sin(thet_tel)*cos(phi_tel) + coord_X_previous;
    vect.y_par = t_var*sin(thet_tel)*sin(phi_tel) + coord_Y_previous;
    vect.z_par = t_var*cos(thet_tel) + coord_Z_previous;

    return vect;
}


angles_for_peel Photon_packet::get_angles_for_peel()
{
    angles_for_peel angles_for_p;
    
    angles_for_p.Thet_peeling = Thet_peel;
    angles_for_p.Thet_teles = thet_tel;
    angles_for_p.Phi_teles = phi_tel;
    angles_for_p.Thet_incid = THET_1;
    angles_for_p.Phi_incid = PHI_1;

    return angles_for_p;
}


void Photon_packet::set_stokes_q_tel(double q_tel)
{
    stokes_q_tel = q_tel; 
}

    
void Photon_packet::set_stokes_u_tel(double u_tel)
{
    stokes_u_tel = u_tel;
}


void Photon_packet::set_stokes_v_tel(double v_tel)
{
    stokes_v_tel = v_tel;
}


void Photon_packet::set_cos_2psi_tel(float cos_2psi_t)
{
    cos_2psi_tel = cos_2psi_t;
}


void Photon_packet::set_sin_2psi_tel(float sin_2psi_t)
{
    sin_2psi_tel = sin_2psi_t;
}