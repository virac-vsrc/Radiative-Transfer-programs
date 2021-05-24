/**
 * @file      photon_packet.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class that simulates the photon packet and its travelling 
 *            in dust cloud.
 *
 * @section   DESCRIPTION
 *
 *            This class essential to calculate optical and geometrical way
 *            for photon packet after each scattering. When optical way is 
 *            generated the integral for passed optical way is calculated 
 *            using trapezoid rule, if integral is bigger than generated 
 *            number, calculations are stopped and photon packet scattering
 *            occurs. The integral calculations using trapezoidal rule were
 *            compared with analytical solution and the precision of 
 *            calculations is high, the average mistake is less than 10e-5.
 *            The photon packet travels inside cloud, if the photon packet 
 *            is out of cloud the calculations should be stopped.
 *
 */

#ifndef PHOTON_PACKET_H_INCLUDED
#define PHOTON_PACKET_H_INCLUDED

#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h> 
#include <time.h>
#include "photon_source.h"
#include "sub_Volume.h"
#include "dust_cloud.h"
#include "ccd_matr.h"
#include "telescope.h"

#define PLANCK_CONST 6.626176e-34  ///< Planck constant
#define c_0 299792458              ///< speed of light in Vacuum



struct photon_vect    ///< structure to calculate photon packet movement
{
    double x_par;     ///< x parameter for photon packet movement
    double y_par;     ///< y parameter for photon packet movement
    double z_par;     ///< z parameter for photon packet movement
};



struct angles_for_peel   ///< structure for class Angles_montecarlo, for peeling off method
{
    double Thet_peeling; ///< Scattering angle in telescope direction
    float Thet_teles;    ///< Polar angle in telescope direction
    float Phi_teles;     ///< Azimuth angle in telescope direction
    float Thet_incid;    ///< Polar angle of incident photon package
    float Phi_incid;     ///< Azimuth angle of incident photon package
};


class Photon_packet
{
    public:

    Photon_packet();


	/** 
    *   @brief   Sets wave length for photon packet.
    *
	*   @param   double length - the length for photon packet
	*
    *   @return  void
    */  
	void set_wave_length(double length);
	
	
	/** 
    *   @brief   Gets wave length of photon packet.
    *  
    *   @return  double - wave length of photon packet
    */  
    double get_wave_length();

	
	/** 
    *   @brief   Sets Stokes vector I parameter.
    *
    *   @param   double I_par - Stokes vector I parameter
    *	
    *   @return  void
    */  
    void set_stokes_I_param(double I_par);
    
	
	/** 
    *   @brief   Gets Stokes vector I parameter.   
    *	
    *   @return  double - Stokes vector I parameter
    */  
	double get_stokes_I_param();
	
	
	/** 
    *   @brief   Sets Stokes vector Q parameter divided by I parameter.
    *
    *   @param   double q_par - Stokes Q parameter divided by I parameter
    *	
    *   @return  void
    */  
    void set_stokes_q_param(double q_par);
	
	
	/** 
    *   @brief   Gets Stokes vector Q parameter divided by I parameter.
    *
    *   @return  double - Q parameter divided by I parameter.
    */  
    double get_stokes_q_param();
	
	
	/** 
    *   @brief   Sets Stokes vector U parameter divided by I parameter.
    *
    *   @param   double u_par - Stokes U parameter divided by I parameter
    *	
    *   @return  void 
    */  
    void set_stokes_u_param(double u_par);
	
	
	/** 
    *   @brief   Gets Stokes vector U parameter divided by I parameter.
    *	
    *   @return  double - Stokes U parameter divided by I parameter
    */  
    double get_stokes_u_param();
	
	
	/** 
    *   @brief   Sets Stokes vector V parameter divided by I parameter.
    *
    *   @param   double v_par - Stokes V parameter divided by I parameter
    *	
    *   @return  void
    */  
    void set_stokes_v_param(double v_par);
	
	
	/** 
    *   @brief   Gets Stokes vector V parameter divided by I parameter.
    *
    *   @return  double - Stokes V parameter divided by I parameter
    */  
    double get_stokes_v_param();

	
	/** 
    *   @brief   Sets optical way for photon packet.
    *
    *   @param   float opt_leng - length of optical way.
    *	
    *   @return  void    
    */  
    void set_optical_way(float opt_leng);
	
	
	/** 
    *   @brief   Gets optical way for photon packet.
    *	
    *   @return  float - length of optical way
    */  
    float get_optical_way();

	
	/** 
    *   @brief   Sets variable that describes photon packet condition.
	*            0 - after emission, 1 - in cloud, 3 - out of cloud, 
    *            4 - energy is less then one photon energy
    *
    *   @param   short location - variable that describes photon 
	*            packet condition
    *	
    *   @return  void
    */  
    void set_photon_location(short location);
	
	
	/** 
    *   @brief   Gets variable that describes photon packet condition.
	*            0 - after emission, 1 - in cloud, 3 - out of cloud, 
    *            4 - photon packet weight is less then one photon 
    *	
    *   @return  short - variable that describes photon packet condition.
    */  
    short get_photon_location();

	
	/** 
    *   @brief   Sets x coordinate of photon packet current position.
    *
    *   @param   float x_cor - x coordinate of photon packet current
	*            position
    *	
    *   @return  void
    */  
    void set_coord_X_current(float x_cor);
	
	
	/** 
    *   @brief   Gets x coordinate of photon packet current position.
    *	
    *   @return  float - x coordinate of photon packet current
	*            position
    */  
    float get_coord_X_current();
	
	
	/** 
    *   @brief   Sets y coordinate of photon packet current position.
    *
    *   @param   float y_cor - y coordinate of photon packet current
	*            position
    *	
    *   @return  void
    */  
    void set_coord_Y_current(float y_cor);
    
	
	/** 
    *   @brief   Gets y coordinate of photon packet current position.
    *	
    *   @return  float - y coordinate of photon packet current
	*            position
    */  
	float get_coord_Y_current();
	
	
	/** 
    *   @brief   Sets z coordinate of photon packet current position.
    *
    *   @param   float z_cor - z coordinate of photon packet current
	*            position
    *	
    *   @return  void
    */  
    void set_coord_Z_current(float z_cor);
	
	
	/** 
    *   @brief   Gets z coordinate of photon packet current position.
    *	
    *   @return  float - z coordinate of photon packet current
	*            position
    */  
    float get_coord_Z_current();

	
	/** 
    *   @brief   Sets x coordinate of photon packet previous position.
    *
    *   @param   float x_coor - x coordinate of photon packet previous
	*            position
    *	
    *   @return  void
    */  
    void set_coord_X_previous(float x_coor);
	
	
	/** 
    *   @brief   Gets x coordinate of photon packet previous position.
    *
    *   @return  float - x coordinate of photon packet previous position.
    */  
    float get_coord_X_previous();
	
	
	/** 
    *   @brief   Sets y coordinate of photon packet previous position.
    *
    *   @param   float y_coor - y coordinate of photon packet previous
	*            position
    *	
    *   @return  void
    */  
    void set_coord_Y_previous(float y_coor);
	
	
	/** 
    *   @brief   Gets y coordinate of photon packet previous position.
    *
    *   @return  float - y coordinate of photon packet previous position.
    */  
    float get_coord_Y_previous();
	
	
	/** 
    *   @brief   Sets z coordinate of photon packet previous position.
    *
    *   @param   float z_coor - z coordinate of photon packet previous
	*            position
    *	
    *   @return  void
    */  
    void set_coord_Z_previous(float z_coor);
	
	
	/** 
    *   @brief   Gets z coordinate of photon packet previous position.
    *
    *   @return  float - z coordinate of photon packet previous position.
    */  
    float get_coord_Z_previous();

	
	/** 
    *   @brief   Sets weight of photon packet.
    *
    *   @param   double phot_weight - weight of photon packet
    *	
    *   @return  void
    */  
    void set_photon_weight(double phot_weight);
	
	
	/** 
    *   @brief   Gets weight of photon packet.
    *	
    *   @return  double - weight of photon packet
    */  
    double get_photon_weight();
	
	
	/** 
    *   @brief   Sets photon packet energy.
    *
    *   @param   double photon_energy - photon packet energy
    *	
    *   @return  void
    */  
    void set_photon_energy(double phot_energ);
    
	
	/** 
    *   @brief   Gets photon packet energy.
    *	
    *   @return  double - photon packet energy
    */  
	double get_photon_energy();

	
    //void calc_phot_energ_0();
	
	
	/** 
    *   @brief   Photon packet scattering. (Reducing photon packet weight)
    *
    *   @param   double albedo - albedo of dust cloud particles
    *	
    *   @return  void
    */  
    void sca_phot_pack(double albedo);

	
	/** 
    *   @brief   Sets azimuth angle of photon packet.
    *
    *   @param   double s_PHI_1 - azimuth angle
    *	
    *   @return  void
    */  
    void set_PHI_1(double s_PHI_1); 
	
	
	/** 
    *   @brief   Gets azimuth angle of photon packet.
    *	
    *   @return  double - azimuth angle
    */  
    double get_PHI_1();
	
	
	/** 
    *   @brief   Sets polar angle of photon packet.
    *
    *   @param   double s_THET_1 - polar angle
    *	
    *   @return  void
    */  
    void set_THET_1(double s_THET_1);
	
	
	/** 
    *   @brief   Gets polar angle of photon packet.
    *	
    *   @return  double - polar angle
    */  
    double get_THET_1();

	
	/** 
    *   @brief   Initializes photon packet with parameters after photon
	*            packet emission.
    *
    *   @param   emit_phot_param pho_struct - structure that contains 
	*            polar and azimuth angles, Stokes parameters, photon
	*            weight, wave length and optical way
    *	
    *   @return  void
    */  
    void initialize_packet(emit_phot_param pho_struct);

    
	/** 
    *   @brief   Calculates the maximum optical way of photon packet 
	*            after emission.    
    *
	*   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  float - maximum optical way after emission. 
    */  
    float max_optical_way_calc(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z]);
    
	
	/** 
    *   @brief   Calculates three dimensional parametric equation of 
	*            line. (using polar and azimuth angles)
    *
    *   @param   float t_var - the parameter for equation
    *	
    *   @return  photon_vect - the coordinates after equation 
	*            calculation.
    */  
    photon_vect vector_calc(float t_var);

    
	/** 
    *   @brief   Calculates geometrical way for photon packet.
    *
    *   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  double - geometrical way
    */  
    double calc_geometrical_way(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z]);

    
    /** 
    *   @brief   Calculating optical way integral using trapezoidal rule
    *
    *   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  double - optical way integral
    */  
	double add_optical_integ_trapez(Sub_Volume *pointe, double optical_integ);


    /** 
    *   @brief   Generates pseido-random number optical way for first scattering.
    *	
    *   @return  float - optical way
    */ 
    float gener_first_optical_way();


	/** 
    *   @brief   Generates pseido-random number optical way.
    *	
    *   @return  float - optical way
    */  
    float gener_optical_way();

    
	/** 
    *   @brief   First scattering of photon packet.
    *
    *   @param   double albedo - the albedo of dust particles.
    *	
    *   @return  void
    */  
    void first_scatter_phot_pack(double albedo);

	
	/** 
    *   @brief   Set cloud type. 1 - the cloud is sphere.
    *
    *   @param   short cloud_typ - the type of cloud
    *	
    *   @return  void
    */  
    void set_cloud_type(short cloud_typ);

	
	/** 
    *   @brief   Set concentration type. 1 - concentration = const,  
	*            2 - concentration = concentration / radius
    *            3 - concentration = concentration / radius^2
    *
    *   @param   short concent_type - concentration type.
    *	
    *   @return  void
    */  
    void set_concentration_type(short concent_type);

	
	/** 
    *   @brief   Calculates optical way using Simpson's method
    *
    *   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  void
    */  
    void check_opt_way_simps(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z]);

	
	/** 
    *   @brief   Calculates optical way using analytical solution.
    *
    *   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  void
    */  
    void check_opt_way_analytical(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z]);

	
	/** 
    *   @brief   Calculates polar and azimuth angles from cloud center.
    *
    *   @param   double x_1 - x coordinate of photon packet
	*   @param   double y_1 - y coordinate of photon packet
	*   @param   double z_1 - z coordinate of photon packet
    *	
    *   @return  void
    */  
    void calc_PHI_THET_central(double x_1, double y_1, double z_1);

	
	/** 
    *   @brief   Sets pointer for virtual CCD matrix.
    *
    *   @param   CCD_matr *matr_p - the pointer for virtual CCD matrix
    *	
    *   @return  void
    */  
    void set_CCD_matr_point(CCD_matr *matr_p);


    /** 
    *   @brief   Calculates coordinates in cloud assymetrical_1, esential to check
    *            if photon package is inside cloud.
    *	
    *   @return  coordinates of vector
    */  
    photon_vect vector_assymetrical_1();


    /** 
    *   @brief   Sets dust cloud radiuses (inner and outer). 
    *
    *   @param   short outer_Rad - the outer radius of dust cloud
    *   @param   short inner_Rad - the inner radius of dust cloud
    *	
    *   @return  void
    */  
    void set_outer_inner_radius(short outer_Rad, short inner_Rad);


    /** 
    *   @brief   Calculates polar and azimuth angles for outer surface of assymetrical_1 cloud.
    *
    *   @param   double x_1 - x coordinate of photon packet
	*   @param   double y_1 - y coordinate of photon packet
	*   @param   double z_1 - z coordinate of photon packet
    *	
    *   @return  void
    */  
    void calc_PHI_THET_assymetrical_1(double x_1, double y_1, double z_1);


    /** 
    *   @brief   Calculates polar and azimuth angles for outer surface of assymetrical_1 cloud.
    *
    *   @param   double precision - how accurate are calculated coordinates
    *   @param   double x_1 - x coordinate of photon packet
	*   @param   double y_1 - y coordinate of photon packet
	*   @param   double z_1 - z coordinate of photon packet
    *	
    *   @return  void
    */  
    void calc_assym_1_rad_angles(double precision, double x_1, double y_1, double z_1);


    /** 
    *   @brief   Initialize star_coordinates **struct_star_p pointer, that represents the
    *            the stars parameters, which are inside dust cloud.
    *
    *   @param   star_coordinates *star_p[star_quantity] - pointer to pointer array of structures
    *            that represents the stars, which are inside dust cloud.
    *	
    *   @return  void
    */  
    void initialize_struct_star_p(star_coordinates *star_p[star_quantity]);


    /**
     *  @brief   Initialize parameter star_numb, which contais the number of photon emitting
     *           star. It is essential to know which star (in current iteration) is emitting
     *           photons packages. 
     * 
     *  @param   short star_n - the number of emitting star.
     * 
     *  @return  void
     */
    void set_star_numb(short star_n);


    /**
     *  @brief   Sets telescope parameters - solid angle, x, y, z, Cartesian coordinates. 
     * 
     *  @param   tel_param tel_p - structure that contains telescope solid angle, x, y, z 
     *                             Cartesian coordinates.
     * 
     *  @return  void
     */
    void set_tel_parameters(tel_param tel_p);


    /**
     *  @brief   Calculates Thet_peel - scattering angle in telescope direction. 
     * 
     *  @return  void
     */
    void calc_Thet_peel();


    /**
     *  @brief   Sets F11 value for angle Thet_peel. 
     * 
     *  @param   float F11_val - value of F11 parameter for angle Thet_peel
     * 
     *  @return  void
     */
    void set_F11_value(float F11_val);


    /**
     *  @brief   Sets F12 value for angle Thet_peel. 
     * 
     *  @param   float F12_val - value of F12 parameter for angle Thet_peel
     * 
     *  @return  void
     */
    void set_F12_value(float F12_val);


    /**
     *  @brief   Calculates probability that photon package is in telescope. 
     * 
     *  @return  void
     */
    void calc_Prob_phot_tel();


    /** 
    *   @brief   Calculates the optical way from scattering point till end
    *            of cloud in telescope direction
    *
	*   @param   Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z] - 
	*            the pointer to three dimensional array of dust cloud.
    *	
    *   @return  float - optical way from scattering point till end of cloud
    *                    in telescope direction. 
    */  
    float optical_way_peel_calc(Sub_Volume *arr[ARRAY_X][ARRAY_Y][ARRAY_Z]);


    /** 
    *   @brief   Calculates three dimensional parametric equation of line.
	*            (using polar and azimuth angles for peeling off method)
    *
    *   @param   float t_var - the parameter for equation
    *	
    *   @return  photon_vect - the coordinates after equation 
	*            calculation.
    */  
    photon_vect vector_calc_peel( float t_var);


    /** 
    *   @brief   Calculates spherical coordinates in direction from scattering
	*            point to telescope.
    *	
    *   @return  void 
    */ 
    void calc_spher_coor_cur_tel();


    /** 
    *   @brief   Gets angles for peeling off method.
    *	
    *   @return  angles_for_peel - angles for peeling off method:
    *            Scattering angle in telescope direction;
    *            Polar angle in telescope direction
    *            Azimuth angle in telescope direction;
    *            Polar angle of incident photon package;
    *            Azimuth angle of incident photon package.
    */ 
    angles_for_peel get_angles_for_peel();


    /** 
    *   @brief   Sets normalized Stokes vector Q parameter
    *	
    *   @param   double q_tel - normalized Stokes vector Q parameter
    * 
    *   @return  void
    */ 
    void set_stokes_q_tel(double q_tel);


    /** 
    *   @brief   Sets normalized Stokes vector U parameter
    *	
    *   @param   double u_tel - normalized Stokes vector U parameter
    * 
    *   @return  void
    */ 
    void set_stokes_u_tel(double u_tel);


    /** 
    *   @brief   Sets normalized Stokes vector V parameter
    *	
    *   @param   double v_tel - normalized Stokes vector V parameter
    * 
    *   @return  void
    */ 
    void set_stokes_v_tel(double v_tel);


    /** 
    *   @brief   Sets angle cos(2*psi_tel)
    *	
    *   @param   float cos_2psi_t - angle cos(2*psi_tel)
    * 
    *   @return  void
    */ 
    void set_cos_2psi_tel(float cos_2psi_t);


    /** 
    *   @brief   Sets angle sin(2*psi_tel)
    *	
    *   @param   float sin_2psi_t - angle sin(2*psi_tel)
    * 
    *   @return  void
    */ 
    void set_sin_2psi_tel(float sin_2psi_t);


    private:
    double wave_length;        ///< the length of wave that emits photon source

    double stokes_I_param;     ///<  Stokes parameter I (for incident and scattered light)
    double stokes_q_param;     ///<  Stokes parameter Q divided by I (for incident and scattered light)
    double stokes_u_param;     ///<  Stokes parameter U divided by I (for incident and scattered light)
    double stokes_v_param;     ///<  Stokes parameter U divided by I (for incident and scattered light)

    double PHI_1;             ///< azimuthal incident angle
    double THET_1;            ///< polar incident angle

    float coord_X_current;    ///< current photon_packet x coordinate inside Sub_volume
    float coord_Y_current;    ///< current photon_packet y coordinate inside Sub_volume
    float coord_Z_current;    ///< current photon_packet z coordinate inside Sub_volume
    float coord_X_previous;   ///< previous photon_packet x coordinate inside Sub_volume
    float coord_Y_previous;   ///< previous photon_packet y coordinate inside Sub_volume
    float coord_Z_previous;   ///< previous photon_packet z coordinate inside Sub_volume

    float prev_x;             ///< to correctly calculate optical_way_integ
    float prev_y;             ///< to correctly calculate optical_way_integ
    float prev_z;             ///< to correctly calculate optical_way_integ
        
    float geometrical_way;    ///< geometrical way after emission or scattering
 
    float optical_way_0;      ///< maximal optical way integral from R_in till R_out 
    float optical_way;        ///< length of optical way for particle after scattering or emitting

    short photon_location;    ///< the condition of photon packet
	                          ///< 0 - after emission, 1 - in cloud
                              ///< 3 - out of cloud, 
                              ///< 4 - photon packet weight is less then one photon

    double photon_weight;      ///< the quantity of photons in one packet
    double photon_energy;      ///< the energy of photon packet

    short inner_Radius;       ///< the inner radius of the cloud
    short outer_Radius;       ///< the outer radius of the cloud

    short cloud_type;         ///< the type of cloud, sphere or another figure , 0 -none (error), 1 - sphere

    short concentration_type; ///< the type of dust particle concentration in a cloud, 0 - none (error)
                              ///< 1 - concentration = const,  2 - concentration = concentration / radius
                              ///< 3 - concentration = concentration / radius^2

    int inner_surface;        ///< how many times optical_way_0 was calculated in inner surface
    int outer_surface;        ///< how many times optical_way_0 was calculated in outer surface
    int inside_cloud;         ///< how many times optical_way_0 was calculated inside the cloud

    double THET_0;            ///< polar angle from center of cloud to previous coordinate
    double PHI_0;             ///< azimuthal angle from center of cloud to previous coordinate

    double THET_0_assym_1;    ///< polar angle for outer surface of assymetrical_1 cloud
    double PHI_0_assym_1;     ///< azimuth angle for outer surface of assymetrical_1 cloud

    double x_par_assym;       ///< recalculated x_par for assym_1 for THET_0_assym_1 and PHI_0_assym_1 calculation
    double y_par_assym;       ///< recalculated y_par for assym_1 for THET_0_assym_1 and PHI_0_assym_1 calculation
    double z_par_assym;       ///< recalculated z_par for assym_1 for THET_0_assym_1 and PHI_0_assym_1 calculation
    double calcul_rad;        ///< recalculated rad for assym_1;

    double vector_assym_1;    ///< recalculated distance of assymetrical_1 surface (from center)

    CCD_matr *matr_point;     ///< pointer for virtual CCD matrix

    star_coordinates **struct_star_p; ///< array of structures, that represents different stars inside dust cloud
    short star_numb;          ///< the number of the photon emitting star
    short current_surface;    ///< which one of cloud surfaces is used to calculate photon travelling


    // ---- For peeling off method ---- //
    double tel_solid_angle;   ///< Solid angle of telescope in steradians
    double tel_X_coord;       ///< x coordinate of telescope Cartesian coordinate system
    double tel_Y_coord;       ///< y coordinate of telescope Cartesian coordinate system
    double tel_Z_coord;       ///< z coordinate of telescope Cartesian coordinate system

    float optical_way_peel;   ///< optical way integral from scattering point till cloud outer surface 
                              ///< in telescope direction

    double Thet_peel;          ///< Scattering angle in telescope direction
    float F11_value;          ///< Set F11 value for value of Thet_peel
    float F12_value;          ///< set F12 value for value of Thet_peel
    
    long double stokes_q_tel;      ///< Stokes normalized q parameters that is observed by telescope
    long double stokes_u_tel;      ///< Stokes normalized u parameters that is observed by telescope
    long double stokes_v_tel;      ///< Stokes normalized v parameters that is observed by telescope

    float cos_2psi_tel;       ///< cos(2*psi_tel) calculated by angles_montecarlo (for Prob_phot_tel)
    float sin_2psi_tel;       ///< sin(2*psi_tel) calculated by angles_montecarlo (for Prob_phot_tel)

    long double Prob_phot_tel;     ///< calculates probability that photon package is in telescope

    float rad_tel;            ///< radius for optical_way_peel calculations (from scattering point to telescope)
    float phi_tel;            ///< azimuth angle for optical_way_peel calculations (from scattering point to telescope)
    float thet_tel;           ///< polar angle for optical_way_peel calculations   (from scattering point to telescope)

    long double phot_weight_tel;   ///< photon weight accumulated by telescope matrix
    //----------------------------------//

    //--- for changing telescope position ---//
    double polar_tel_pos;   ///< Telescope position polar angle
    double azimuth_tel_pos; ///< Telescope position azimuthal angle

    double x_tel_vect;      ///< x vector of telescope's XY plane position
    double y_tel_vect;      ///< y vector of telescope's XY plane position
    double z_tel_vect;      ///< z vector of telescope's XY plane position

    double distance_tel_0;  ///< distance between cloud system and telescope system
    double distance_tel;    ///< distance after telescopes XY plane rotations (rotations was made by telescope class)
    //----------------------------------//

};

#endif