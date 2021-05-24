/**
 * @file      sub_Volume.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Sub_Volume class functions.
 */


#include "sub_Volume.h"


Sub_Volume::Sub_Volume()
{

    /*/ for homogen cloud
    scattering_cross_sec = 0.12067;
    extinction_cross_sec = 0.17239;
    //part_concentration = 4; // 20 for homogen cloud
    part_concentration = 2; // 10 for homogen cloud
    //part_concentration = 1; // 5 for homogen cloud
    /*/


    /*/ for  concentration / radius
    scattering_cross_sec = 0.12854;
    extinction_cross_sec = 0.18363;
    //part_concentration = 80; // 20 optical way
    part_concentration = 40; // 10 optical way
    //part_concentration = 20; // 5 optical way
    /*/

    // for   concentration / radius^2
    //scattering_cross_sec = 0.1343967;
    //extinction_cross_sec = 0.1919953;
    //part_concentration = 1400; // 20 optical way
    //part_concentration = 700; // 10 optical way
    //part_concentration = 350; // 5 optical way
    //

    // for homogen cloud
    //scattering_cross_sec = 1.7376;
    //extinction_cross_sec = 1.871;
    //part_concentration = 1.843;   // 100 optical way
    //part_concentration = 0.1843;   // 10 optical way
    //part_concentration = 0.01843;  // 1 optical way
    //part_concentration = 0.001843;   // 0.1 optical way

    // for homogen cloud for analytical program
    //scattering_cross_sec = 0.9355; // albedo = 0.5
    scattering_cross_sec = 1.7376; // albedo = 0.9287
    //scattering_cross_sec = 1.871;  // albedo = 1
    extinction_cross_sec = 1.871;
    //part_concentration = 0.098455;  // 7 optical way
    //part_concentration = 0.08439;   // 6 optical way
    //part_concentration = 0.070325;  // 5 optical way
    part_concentration = 0.014065;  // 1 optical way
    //part_concentration = 0.0014065;  // 0.1 optical way

}

Sub_Volume::~Sub_Volume()
{}




void Sub_Volume::set_extinction(float extinction)
{
    extinction_cross_sec = extinction;
}

float Sub_Volume::get_extinction()
{
    return extinction_cross_sec;
}

void Sub_Volume::set_scattering(float scattering)
{
    scattering_cross_sec = scattering;
}

float Sub_Volume::get_scattering()
{
    return scattering_cross_sec;
}


void Sub_Volume::set_concentration(float number)
{
    part_concentration = number;
}


float Sub_Volume::Sub_Volume::get_concentration()
{
    return part_concentration;
}


void Sub_Volume::set_cloud_surface(short cloud_surf)
{
    cloud_surface = cloud_surf;
}


short Sub_Volume::get_cloud_surface()
{
    return cloud_surface;
}
