/**
 * @file      dust_cloud.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class to create the shape of dust cloud using 
 *            3-dimensional array.
 *	       
 *
 * @section   DESCRIPTION
 *
 *            The class is used to create the dust cloud in 3-dimensional
 *            array with dimensions [ARRAY_X][ARRAY_Y][ARRAY_Z] and type 
 *            short. When minimum and maximum azimuth angle phi, polar 
 *            angle thet also outer and inner radiuses are set, the inner
 *            and outer surfaces are created. The volume of dust cloud 
 *            that is between inner and outer surfaces is filled . This
 *            array is essential to create a 3-dimensional array of 
 *            pointers, which are pointing to dust cloud Sub_Volume 
 *            objects.
 */
 
#ifndef DUST_CLOUD_H_INCLUDED
#define DUST_CLOUD_H_INCLUDED

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265359   ///< PI number

#define ARRAY_X 80         ///< X dimensions for 3D array
#define ARRAY_Y 80         ///< Y dimensions for 3D array
#define ARRAY_Z 80         ///< Z dimensions for 3D array

#define star_quantity 10   ///< array size, how many stars could be used,
                           ///< should be more or equeal to star quantity

using namespace std;


//------------------------How to use class Dust_cloud------------------------------//
//                                                                                 //
//      1. Create the object, use "1" to set the cloud as a sphere,                //
//         use "2" to set the cloud as assymetrical_1.                             //
//                                                                                 //
//      2. Set minimum and maximum phi-azimuth and thet-polar angles for           //
//         surface parametric functions, using set_min_max_phi_thet().             //
//                                                                                 //
//      3. Set inner and outer radius for surface parametric functions, using      //
//         set_big_and_small_rad().                                                //
//                                                                                 //
//      4. Create the shape of cloud using function make_cloud_surface().          //
//                                                                                 //
//      5. Cut off spheres from the dust cloud (inside those spheres dust do not   //
//         exist, only star with radius "star_coordinates.radius"), using          //
//         function cut_off_sphere().                                              //
//                                                                                 //
//      6. Fill dust cloud fill_cloud().                                           //
//                                                                                 //
//      7. It is possible to read the values of 3d array using function            //
//         get_array_value(), without changing any element of an array.            //
//                                                                                 //
//---------------------------------------------------------------------------------//


struct star_coordinates
{
    float radius;       // radius of star
    float cloud_radius; // radius around star, there are no cosmic dust
    float x_pos;        // x - position of star
    float y_pos;        // y - positien of star
    float z_pos;        // z - position of star
};


class Dust_cloud
{
    public:
    Dust_cloud(short cloud_Type);
    ~Dust_cloud();

	
	/** 
    *   @brief   Parametric function which calculates x coordinates
    *            for sphere surface.	
    *  
    *   @param   double phi - azimuth angle for parametric function
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - x coordinate for sphere surface
    */  
    double sphere_X_func(double phi, double thet);
	
	
	/** 
    *   @brief   Parametric function which calculates y coordinates
    *            for sphere surface.	
    *  
    *   @param   double phi - azimuth angle for parametric function
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - y coordinate for sphere surface
    */  
    double sphere_Y_func(double phi, double thet);
	
	
	/** 
    *   @brief   Parametric function which calculates z coordinates
    *            for sphere surface.	
    *  
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - z coordinate for sphere surface
    */  
    double sphere_Z_func(double thet);


    /** 
    *   @brief   Parametric function which calculates x coordinates
    *            for assymetrical_1 surface.	
    *  
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - x coordinate for assymetrical_1 surface
    */  
    double assymetrical_1_X_func(double thet);
	
	
	/** 
    *   @brief   Parametric function which calculates y coordinates
    *            for assymetrical_1 surface.	
    *  
    *   @param   double phi - azimuth angle for parametric function
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - y coordinate for assymetrical_1 surface
    */  
    double assymetrical_1_Y_func(double phi, double thet);
	
	
	/** 
    *   @brief   Parametric function which calculates z coordinates
    *            for assymetrical_1 surface.	
    *
    *   @param   double phi - azimuth angle for parametric function  
	*   @param   double thet - polar angle for parametric funcion
    *
	*   @return  double - z coordinate for assymetrical_1 surface
    */  
    double assymetrical_1_Z_func(double phi, double thet);

    
	/** 
    *   @brief   Sets minimum and maximum azimuth and polar angle values
    *            in radians.	
    *  
    *   @param   double min_Phi - minimum value of azimuth angle
	*   @param   double min_Thet - minimum value of polar angle
	*   @param   double max_Phi - maximum value of azimuth angle
	*   @param   double max_Thet - maximum value of polar angle
    *
	*   @return  void
    */  
    void set_min_max_phi_thet(double min_Phi, double min_Thet, double max_Phi, double max_Thet);
    
    
	/** 
    *   @brief   Sets inner and outer radiuses for dust cloud.
    *  
    *   @param   double out_Rad - outer radius of dust cloud
    *
	*   @return  void
    */  
    void set_outer_rad(double out_Rad);

    
	/** 
    *   @brief   Gets 3-dimensional array element value. 1 - outer surface
    *            2 - volume between surfaces, 3 - inner surface
    *  
    *   @param   short i - first coordinate of element in array
	*   @param   short j - second coordinate of element in array
	*   @param   short k - third coordinate of element in array
    *
	*   @return  short - value of 3-dimensional array element
    */  
    short get_array_value(short i, short j, short k);
    
    
	/** 
    *   @brief   Makes cloud outer surface.
    *
	*   @return  void
    */  
    void make_cloud_surface();
    
	
	/** 
    *   @brief   Makes cloud inner surface, by cutting off spheres from 
	*            dust cloud (inside those spheres no dust particles only stars
    *            with a radius "star_coordinates.radius").
    *
    *   @param   star_coordinates *st_coord[10] - array of stars, each element has:
    *            radius of sphere around star, and it's x,y,z position.
    * 
	*   @return  void
    */  
    void cut_off_sphere(star_coordinates *st_coord[star_quantity]);

    
	/** 
    *   @brief   Fils the volume of dust cloud between outer and inner
	*            surfaces.
    *
	*   @return  void
    */  
    void fill_cloud();


    private:

    short three_dim_array[ARRAY_X][ARRAY_Y][ARRAY_Z];  ///< 3-dimensional array for dust cloud    

    short cloud_type;                                  ///< variable to select the shape of dust cloud

    double outer_rad;                                  ///< outer radius of cloud for cloud shape functions

    double min_phi;                                    ///< minimum azimuth angle Phi for cloud shape functions
    double max_phi;                                    ///< maximum azimuth angle Phi for cloud shape functions
    double min_thet;                                   ///< minimum polar angle Thet for cloud shape functions
    double max_thet;                                   ///< maximum polar angle Thet for cloud shape functions

    short cloud_UP;                                    ///< variable to check if above the 3-d array element is cloud surface
    short cloud_DOWN;                                  ///< variable to check if below the 3-d array element is cloud surface
    short cloud_FRONT;                                 ///< variable to check if in front of the 3-d array element is cloud surface
    short cloud_BACK;                                  ///< variable to check if behind the 3-d array element is cloud surface
    short cloud_RIGHT;							       ///< variable to check if on the right of the 3-d array element is cloud surface
    short cloud_LEFT;                                  ///< variable to check if on the left of the 3-d array element is cloud surface
};

#endif