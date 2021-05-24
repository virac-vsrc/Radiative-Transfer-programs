/**
 * @file      telescope.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Telescope class functions.
 */

#include "telescope.h"


Telescope::Telescope()
{
    cloud_radius = -1;      
    dist_cloud_telesc = -1; 
    telescope_rad = -1;     
    telescope_area = -1;    
    solid_angle = -1; 
}


Telescope::~Telescope()
{}


void Telescope::set_cloud_rad_meters(double meters)
{
    cloud_radius = meters;
}


double Telescope::get_cloud_rad_meters()
{
    return cloud_radius;
}


void Telescope::set_cloud_rad_AU(double astron_units)
{
    cloud_radius = astron_units*AU;
}


double Telescope::get_cloud_rad_AU()
{
    return cloud_radius/AU;
}


void Telescope::set_dist_cloud_tel_meters(double meters)
{
    dist_cloud_telesc = meters;
}


double Telescope::get_dist_cloud_tel_meters()
{
    return dist_cloud_telesc;
}


void Telescope::set_dist_cloud_tel_PC(double parsecs)
{
    dist_cloud_telesc = parsecs*PC;
}


double Telescope::get_dist_cloud_tel_PC()
{
    return dist_cloud_telesc/PC;
}


void Telescope::set_telescope_rad_meters(double meters)
{
    telescope_rad = meters;
}


double Telescope::get_telescope_rad_meters()
{
    return telescope_rad;
}


void Telescope::calc_telescope_area()
{
    telescope_area = telescope_rad*telescope_rad*PI;
}


double Telescope::get_telescope_area()
{
    return telescope_area;
}


void Telescope::calc_solid_angle()
{
    solid_angle = telescope_area / (dist_cloud_telesc*dist_cloud_telesc);
}


double Telescope::get_solid_angle()
{
    return solid_angle;
}


void Telescope::set_tel_X_coord(double x)
{
    tel_X_coord = x;
}


double Telescope::get_tel_X_coord()
{
    return tel_X_coord;
}
    

void Telescope::set_tel_Y_coord(double y)
{
    tel_Y_coord = y;
}

    
double Telescope::get_tel_Y_coord()
{
    return tel_Y_coord;
}


void Telescope::set_tel_Z_coord(double z)
{
    tel_Z_coord = z;
}

    
double Telescope::get_tel_Z_coord()
{
    return tel_Z_coord;
}


void Telescope::set_tel_polar_angle(double thet_angle)
{
    if(thet_angle < 0 || thet_angle > PI)
    {
        cout << "error - telescope's position polar angle should be bigger or equal to 0 and less or equal to PI.\n";
    }

    tel_polar_angle = thet_angle;       
}


void Telescope::set_tel_azimuth_angle(double phi_angle)
{
    if(phi_angle < 0 || phi_angle > 2*PI)
    {
        cout << "error - telescope's position azimuthal angle should be bigger or equal to 0 and less than 2*PI.\n";
    }

    tel_azimuth_angle = phi_angle;
}


void Telescope::set_tel_dist_cloud_sys_tel_sys(double distance)
{
    if(distance < 0)
    {
        cout << "error - distance between cloud system and telescope system should be positive. \n";
    }
    else
    {
        dist_cloud_sys_tel_sys = distance + ARRAY_Z; // the start position of telescope should be above cloud
    }
    
}


void Telescope::calc_unity_vect_tel_XY_plane()
{
    x_vect = sin(tel_polar_angle)*cos(tel_azimuth_angle);
    y_vect = sin(tel_polar_angle)*sin(tel_azimuth_angle);
    z_vect = cos(tel_polar_angle);
}


void Telescope::calc_dist_for_tel_XY_plane()
{
    dist_for_calculations = x_vect*tel_X_coord + y_vect*tel_Y_coord + z_vect*tel_Z_coord;
}


void Telescope::calc_parameters_for_tel_param()
{
    // start position of telescope
    double x_0 = ARRAY_X/2;
    double y_0 = ARRAY_Y/2;
    double z_0 = dist_cloud_sys_tel_sys;

    // changing position before polar angle rotation
    double x_1 = x_0 - ARRAY_X/2;
    double y_1 = y_0;
    double z_1 = z_0 - ARRAY_Z/2;

    // polar angle rotation
    double x_2 = x_1*cos(tel_polar_angle) + z_1*sin(tel_polar_angle);
    double y_2 = y_1;
    double z_2 = -x_1*sin(tel_polar_angle) + z_1*cos(tel_polar_angle);

    // changing position after polar angle rotation
    double x_3 = x_2 + ARRAY_X/2;
    double y_3 = y_2;
    double z_3 = z_2 + ARRAY_Z/2;

    // changing position before azimuthal angle rotation
    double x_4 = x_3 - ARRAY_X/2;
    double y_4 = y_3 - ARRAY_Y/2;
    double z_4 = z_3;

    // azimuthal angle rotation
    double x_5 = x_4*cos(tel_azimuth_angle) - y_4*sin(tel_azimuth_angle);
    double y_5 = x_4*sin(tel_azimuth_angle) + y_4*cos(tel_azimuth_angle);
    double z_5 = z_4;

    // changing position after azimuthal angle rotation
    tel_X_coord = x_5 + ARRAY_X/2;
    tel_Y_coord = y_5 + ARRAY_Y/2;
    tel_Z_coord = z_5;

    // calculate unity vector for telescope's XY plane
    calc_unity_vect_tel_XY_plane();

    // calculate distance for telescope's XY plane
    calc_dist_for_tel_XY_plane();
}


tel_param Telescope::get_tel_param()
{
    tel_par.tel_solid_angle = solid_angle;
    tel_par.tel_X_coor = tel_X_coord;
    tel_par.tel_Y_coor = tel_Y_coord;
    tel_par.tel_Z_coor = tel_Z_coord;
    tel_par.teles_polar_angle = tel_polar_angle;
    tel_par.teles_azimuth_angle = tel_azimuth_angle;
    tel_par.distan_cloud_sys_tel_sys = dist_cloud_sys_tel_sys;
    tel_par.tel_dist_for_calculations = dist_for_calculations;
    tel_par.tel_x_vect = x_vect;
    tel_par.tel_y_vect = y_vect;
    tel_par.tel_z_vect = z_vect;

    return tel_par;
}