/**
 * @file      sub_Volume.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class contains optical parameters of dust cloud 
 *            subvolume.
 *
 * @section   DESCRIPTION
 *
 *            The dust cloud is splitted into elements of 3-dimensional
 *            array. Each array element is the subvolume of dust cloud.
 *            Subvolume contains average scattering cross section, 
 *            average extinction cross section, particle concentration, 
 *            parameter that describes subvolume as part of dust cloud 
 *            (outer surface, inner surface or not surface).
 *            
 */

#ifndef SUB_VOLUME_H_INCLUDED
#define SUB_VOLUME_H_INCLUDED

#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;


class Sub_Volume
{
public:
    Sub_Volume(); // creates the object atomaticly setting object parameters
    ~Sub_Volume();

	
	/** 
    *   @brief   Sets average extinction cross section for dust cloud
    *	         subvolume.   
    *  
    *   @param   float extinction - average extinction cross section
	*
	*   @return  void
    */  
    void set_extinction(float extinction);
	
	
	/** 
    *   @brief   Gets average extinction cross section of dust cloud
    *            subvolume.	
    *
	*   @return  float - average extinction cross section
    */  
    float get_extinction();
	
	
	/** 
    *   @brief   Sets average scattering cross section for dust cloud
    *	         subvolume.   
    *  
    *   @param   float scattering - average scattering cross section
	*
	*   @return  void
    */  
    void set_scattering(float scattering);
	
	
	/** 
    *   @brief   Gets average scattering cross section of dust cloud
    *            subvolume.	
    *
	*   @return  float - average scattering cross section
    */  
    float get_scattering();                

	
	/** 
    *   @brief   Sets the number of particles in cm^3 for subvolume.
    *  
    *   @param   float number - number of particles in cm^3
	*
	*   @return  void
    */  
    void set_concentration(float number);
	
	
	/** 
    *   @brief   Gets the number of particles in cm^3 of subvolume.
	*
	*   @return  float - number of particles in cm^3
    */  
    float get_concentration();            

	
	/** 
    *   @brief   Sets the type of subvolume as a part of dust cloud.
	*            1 - outer surface, 2-not surface (volume between surfaces),
    *            3 - inner surface.
    *  
    *   @param   short cloud_surf - type of subvolume as a part of dust
	*            cloud
	*
	*   @return  void
    */  
    void set_cloud_surface(short cloud_surf);
	

	/** 
    *   @brief   Gets the type of subvolume as a part of dust cloud.
	*            1 - outer surface, 2-not surface (volume between surfaces),
    *            3 - inner surface.
	*
	*   @return  short - type of subvolume as a part of dust cloud
    */  
    short get_cloud_surface();


private:
    float extinction_cross_sec;   ///< extinction cross section
    float scattering_cross_sec;   ///< scattering cross section

    float part_concentration;     ///< the number of particle in cm^3 

    short cloud_surface;          ///< type of subvolume as a part of dust cloud
                                  ///< 1 - outer surface, 2 - not surface, 3 - inner surface
};

#endif // SUB_VOLUME_H_INCLUDED