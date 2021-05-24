/**
 * @file      telescope.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class contains parameters of telescope and 
 *            astronomical constants.
 *
 * @section   DESCRIPTION
 *
 *            The class contains physical dimensions of the
 *            cloud, distance between cloud and telescope,
 *            radius of telescope. Calculates solid angle of
 *            telescope.
 *            
 */

#ifndef TELESCOPE_H_INCLUDED
#define TELESCOPE_H_INCLUDED

#include <iostream>
#include <cmath>
#include "dust_cloud.h"

#define PC 3.0856775815e+16      ///< 1 Parsec in meters
#define AU 1.495978707e+11       ///< 1 Astronomical Unit in meters


struct tel_param             ///< structure of telescope parameters
{
    double tel_solid_angle;  ///< solid angle of telescope

    double tel_X_coor;       ///< x Cartesian coordinate of telescope
    double tel_Y_coor;       ///< y Cartesian coordinate of telescope
    double tel_Z_coor;       ///< z Cartesian coordinate of telescope

    double teles_polar_angle;        ///< Telescope position polar angle
    double teles_azimuth_angle;      ///< Telescope position azimuthal angle
    double distan_cloud_sys_tel_sys; ///< distance between cloud system and telescope system

    double tel_dist_for_calculations;  ///< distance after telescopes XY plane rotation

    double tel_x_vect; ///< x vector of telescope's XY plane position 
    double tel_y_vect; ///< y vector of telescope's XY plane position 
    double tel_z_vect; ///< z vector of telescope's XY plane position 
};


class Telescope
{
public:
    Telescope();
    ~Telescope();


    /** 
    *   @brief   Sets cloud radius using meters.   
    *  
    *   @param   double meters - radius of dust cloud in meters
	*
	*   @return  void
    */  
    void set_cloud_rad_meters(double meters);


    /** 
    *   @brief   Gets cloud radius in meters.   
    *  
    *   @param   none
	*
	*   @return  double - cloud radius in meters
    */  
    double get_cloud_rad_meters();


    /** 
    *   @brief   Sets cloud radius in astronomical units.   
    *  
    *   @param   double astron_units - radius of dust cloud in astronomical units.
	*
	*   @return  void
    */
    void set_cloud_rad_AU(double astron_units);


    /** 
    *   @brief   Gets cloud radius in astronomical units.   
    *  
    *   @param   none
	*
	*   @return  double - radius of dust cloud in astronomical units.
    */
    double get_cloud_rad_AU();


    /** 
    *   @brief   Sets distance between cloud and telescope in meters.   
    *  
    *   @param   double meters - distance between cloud and telescope in meters
	*
	*   @return  void
    */
    void set_dist_cloud_tel_meters(double meters);
    
    
    /** 
    *   @brief   Gets distance between cloud and telescope in meters.   
    *  
    *   @param   none
	*
	*   @return  double - distance between cloud and telescope in meters
    */
    double get_dist_cloud_tel_meters();


    /** 
    *   @brief   Sets distance between cloud and telescope in parsecs.   
    *  
    *   @param   double parsecs - distance between cloud and telescope in parsecs
	*
	*   @return  void
    */
    void set_dist_cloud_tel_PC(double parsecs);
    
    
    /** 
    *   @brief   Gets distance between cloud and telescope in parsecs.   
    *  
    *   @param   none 
	*
	*   @return  double - distance between cloud and telescope in parsecs
    */
    double get_dist_cloud_tel_PC();


    /** 
    *   @brief   Sets radius of telescope in meters.   
    *  
    *   @param   double meters - radius of the telescope in meters
	*
	*   @return  void
    */
    void set_telescope_rad_meters(double meters);


    /** 
    *   @brief   Gets radius of telescope in meters.   
    *  
    *   @param   none
	*
	*   @return  double - radius of the telescope in meters
    */
    double get_telescope_rad_meters();


    /** 
    *   @brief   Calculates telescope area in meters^2.   
    *  
    *   @param   none
	*
	*   @return  void
    */
    void calc_telescope_area();


    /** 
    *   @brief   Gets telescope area in meters^2.   
    *  
    *   @param   none
	*
	*   @return  double - telescope area in meters^2.
    */
    double get_telescope_area();


    /** 
    *   @brief   Calculates solid angle in steradians.   
    *  
    *   @param   none
	*
	*   @return  void
    */
    void calc_solid_angle();


    /** 
    *   @brief   Gets solid angle in steradians.   
    *  
    *   @param   none
	*
	*   @return  double - solid angle in steradians.
    */
    double get_solid_angle();


    /** 
    *   @brief   Sets Cartesian coordinate x of telescope.   
    *  
    *   @param   double x - Cartesian coordinate x of telescope.
	*
	*   @return  void
    */
    void set_tel_X_coord(double x);


    /** 
    *   @brief   Gets Cartesian coordinate x of telescope.   
    *  
    *   @param   none 
	*
	*   @return  double - Cartesian coordinate x of telescope.
    */
    double get_tel_X_coord();
    

    /** 
    *   @brief   Sets Cartesian coordinate y of telescope.   
    *  
    *   @param   double y - Cartesian coordinate y of telescope.
	*
	*   @return  void
    */
    void set_tel_Y_coord(double y);


    /** 
    *   @brief   Gets Cartesian coordinate y of telescope.   
    *  
    *   @param   none 
	*
	*   @return  double - Cartesian coordinate y of telescope.
    */
    double get_tel_Y_coord();


    /** 
    *   @brief   Sets Cartesian coordinate z of telescope.   
    *  
    *   @param   double z - Cartesian coordinate z of telescope.
	*
	*   @return  void
    */
    void set_tel_Z_coord(double z);


    /** 
    *   @brief   Gets Cartesian coordinate z of telescope.   
    *  
    *   @param   none 
	*
	*   @return  double - Cartesian coordinate z of telescope.
    */
    double get_tel_Z_coord();


    /** 
    *   @brief   Sets the telescope's position polar angle.  
    *  
    *   @param   double thet_angle - telescope's position polar angle
	*
	*   @return  void
    */
    void set_tel_polar_angle(double thet_angle);


    /** 
    *   @brief   Sets the telescope's position azimuthal angle.  
    *  
    *   @param   double phi_angle - telescope's position azimuthal angle
	*
	*   @return  void
    */
    void set_tel_azimuth_angle(double phi_angle);


    /** 
    *   @brief   Sets distance between cloud system and telescope system 
    *  
    *   @param   double distance - distance between cloud system and telescope system.
	*
	*   @return  void
    */
    void set_tel_dist_cloud_sys_tel_sys(double distance);


    /** 
    *   @brief   Calculate unity vector of telescope's XY plane.
    *  
    *   @param   void
	*
	*   @return  void
    */
    void calc_unity_vect_tel_XY_plane();


    /** 
    *   @brief   Calculate distance for telescope's XY plane equation.
    *  
    *   @param   void
	*
	*   @return  void
    */
    void calc_dist_for_tel_XY_plane();


    /** 
    *   @brief   Calculate parameters for tel_param structure.
    *  
    *   @param   void
	*
	*   @return  void
    */
    void calc_parameters_for_tel_param();


    /** 
    *   @brief   Gets parameters of telescope.   
    *  
    *   @param   none 
	*
	*   @return  tel_param - structure of telescope parameters.
    */
    tel_param get_tel_param();


private:
    double cloud_radius;      ///< Physical dimensions of cloud in meters
    double dist_cloud_telesc; ///< Distance between cloud and telescope in meters
    double telescope_rad;     ///< Radius of the telescope in meters
    double telescope_area;    ///< Area of the telescope in meters^2
    double solid_angle;       ///< Solid angle of telescope in steradians
    
    double tel_X_coord;       ///< x coordinate of telescope Cartesian coordinate system
    double tel_Y_coord;       ///< y coordinate of telescope Cartesian coordinate system
    double tel_Z_coord;       ///< z coordinate of telescope Cartesian coordinate system
    
    double tel_polar_angle;        ///< Telescope position polar angle
    double tel_azimuth_angle;      ///< Telescope position azimuthal angle
    double dist_cloud_sys_tel_sys; ///< distance between cloud system and telescope system

    double dist_for_calculations;  ///< distance after telescopes XY plane rotation

    double x_vect; ///< x vector of telescope's XY plane position 
    double y_vect; ///< y vector of telescope's XY plane position 
    double z_vect; ///< z vector of telescope's XY plane position 

    tel_param tel_par;        ///< parameters of telescope
};


#endif // TELESCOPE_H_INCLUDED