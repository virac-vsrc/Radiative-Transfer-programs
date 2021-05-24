/**
 * @file      photon_source.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Photon_source class functions.
 */
 

#include "photon_source.h"

using namespace std;

Photon_source::Photon_source()
{
    wave_length = 0; 
    
    stokes_I_param = 0; 
    stokes_Q_param = 0; 
    stokes_U_param = 0; 
    stokes_V_param = 0; 
    
    PHI_1 = 0;           
    THET_1 = 0;          
     
    photon_weight = 0;

    phot_pack_X = 0;
    phot_pack_Y = 0;
    phot_pack_Z = 0;

    srand(time(NULL)); 
}


void Photon_source::set_wave_length(double length)
{
    wave_length = length;
}


double Photon_source::get_wave_length()
{
    return wave_length;
}


void Photon_source::set_stokes_I_param(double I_par)
{
    stokes_I_param = I_par;
}


double Photon_source::get_stokes_I_param()
{
    return stokes_I_param;
}


void Photon_source::set_stokes_Q_param(double Q_par)
{
    stokes_Q_param = Q_par;
}


double Photon_source::get_stokes_Q_param()
{
    return stokes_Q_param;
}


void Photon_source::set_stokes_U_param(double U_par)
{
    stokes_U_param = U_par;
}


double Photon_source::get_stokes_U_param()
{
    return stokes_U_param;
}


void Photon_source::set_stokes_V_param(double V_par)
{
    stokes_V_param = V_par;
}


double Photon_source::get_stokes_V_param()
{
    return stokes_V_param;
}
   

void Photon_source::set_phot_weight(double phot_weight)
{
    photon_weight = phot_weight;
}


void Photon_source::generate_PHI_THET_surf()
{
    //PHI_surf = (rand() % 360000 ) / 1000.0 / 180.0 * PI; // from 0 to 2*PI radians
    PHI_surf = (rand() % 3600 ) / 10.0 / 180.0 * PI;

    //double cos_value = (rand() % 200001) / 100000.0 - 1.0; // from -1 to +1
    double cos_value = (rand() % 2001) / 1000.0 - 1.0;
    THET_surf = acosf(cos_value); // from 0 to PI radians
}


// output azimuthal, polar angles, stokes vector parameters
emit_phot_param Photon_source::emit_phot()
{
    emit_phot_param output;

    output.s_PHI_1 = PHI_1;
    output.s_THET_1 = THET_1;
    
    output.stokes_I = stokes_I_param;
    output.stokes_Q = stokes_Q_param;
    output.stokes_U = stokes_U_param;
    output.stokes_V = stokes_V_param;

    output.phot_weight = photon_weight;
    output.wave_leng = wave_length;

    output.phot_pack_start_X = phot_pack_X;
    output.phot_pack_start_Y = phot_pack_Y;
    output.phot_pack_start_Z = phot_pack_Z;

    return output;
}


void Photon_source::initialize_photon_source(star_coordinates *star_coord)
{
    star_center_X = star_coord->x_pos;
    star_center_Y = star_coord->y_pos;
    star_center_Z = star_coord->z_pos;

    star_radius = star_coord->radius;
}


void Photon_source::gener_emit_angles_point()
{
    generate_PHI_THET_surf();
    calc_emition_point();
    calc_normal_unit_vect();
    gener_delta_THET_alfa();
    calc_emit_0_vect();
    calc_emit_1_vect();
    calc_emition_angles();
}


void Photon_source::calc_emition_point()
{
    phot_pack_X = star_radius*sin(THET_surf)*cos(PHI_surf) + star_center_X;
    phot_pack_Y = star_radius*sin(THET_surf)*sin(PHI_surf) + star_center_Y;
    phot_pack_Z = star_radius*cos(THET_surf) + star_center_Z;
}


void Photon_source::calc_normal_unit_vect()
{
    un_norm_vect_x = sin(THET_surf)*cos(PHI_surf);
    un_norm_vect_y = sin(THET_surf)*sin(PHI_surf);
    un_norm_vect_z = cos(THET_surf);
}


void Photon_source::gener_delta_THET_alfa()
{
    //float cos_value = (rand() % 100001) / 100000.0; // from 0 to +1
    float cos_value = (rand() % 1001) / 1000.0; // from 0 to +1
    delta_THET = acosf(cos_value); // from 0 to PI/2 radians
    //alfa = (rand() % 360000 ) / 1000.0 / 180.0 * PI; // from 0 to 2*PI radians
    alfa = (rand() % 3600 ) / 10.0 / 180.0 * PI; 
}


void Photon_source::calc_emit_0_vect()
{
    un_em_0_vect_x = sin(THET_surf + delta_THET)*cos(PHI_surf);
    un_em_0_vect_y = sin(THET_surf + delta_THET)*sin(PHI_surf);
    un_em_0_vect_z = cos(THET_surf + delta_THET);
}


void Photon_source::calc_emit_1_vect()
{
    un_em_1_vect_x = un_em_0_vect_x*(cos(alfa)+(1-cos(alfa))*un_norm_vect_x*un_norm_vect_x)
                     + un_em_0_vect_y*((1-cos(alfa))*un_norm_vect_x*un_norm_vect_y - sin(alfa)*un_norm_vect_z)
                     + un_em_0_vect_z*((1-cos(alfa))*un_norm_vect_x*un_norm_vect_z + sin(alfa)*un_norm_vect_y);

    un_em_1_vect_y = un_em_0_vect_x*((1-cos(alfa))*un_norm_vect_x*un_norm_vect_y + sin(alfa)*un_norm_vect_z)
                     + un_em_0_vect_y*(cos(alfa) + (1-cos(alfa))*un_norm_vect_y*un_norm_vect_y)
                     + un_em_0_vect_z*((1-cos(alfa))*un_norm_vect_y*un_norm_vect_z - sin(alfa)*un_norm_vect_x);

    un_em_1_vect_z = un_em_0_vect_x*((1-cos(alfa))*un_norm_vect_x*un_norm_vect_z - sin(alfa)*un_norm_vect_y)
                     + un_em_0_vect_y*((1-cos(alfa))*un_norm_vect_y*un_norm_vect_z + sin(alfa)*un_norm_vect_x)
                     + un_em_0_vect_z*(cos(alfa) + (1-cos(alfa))*un_norm_vect_z*un_norm_vect_z);
}


void Photon_source::calc_emition_angles()
{
    THET_1 = acos(un_em_1_vect_z / sqrt(un_em_1_vect_x*un_em_1_vect_x + un_em_1_vect_y*un_em_1_vect_y
                                        + un_em_1_vect_z*un_em_1_vect_z));

    PHI_1 = atan2(un_em_1_vect_y, un_em_1_vect_x);
}