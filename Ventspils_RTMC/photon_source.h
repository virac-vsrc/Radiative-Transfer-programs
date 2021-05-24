/**
 * @file      photon_source.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class to emit photon packets from photon source. 
 *
 * @section   DESCRIPTION
 *
 *			  The class generates polar and azimuth angle for photon packet. 
 *            It is essential to set the Stokes vector parameters and photon
 *            weight. When all parameters of photon packet are set the photon
 *            packet is emitted form photon source.
 */
 
#ifndef PHOTON_SOURCE_H_INCLUDED
#define PHOTON_SOURCE_H_INCLUDED

#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include "dust_cloud.h"


struct emit_phot_param        ///< structure that contains data of emitted photon packet
{
    double s_PHI_1;           ///< azimuthal angle of emitted photon packet
    double s_THET_1;          ///< polar angle of emitted photon packet

    double stokes_I;          ///< Stokes I parameter of photon packet
    double stokes_Q;          ///< Stokes Q parameter of photon packet
    double stokes_U;          ///< Stokes U parameter of photon packet
    double stokes_V;          ///< Stokes V parameter of photon packet

    double phot_weight;       ///< weight of photon packet
    double wave_leng;         ///< the wave length of photon packet

    double phot_pack_start_X; ///< the starting X coordinate of photon package
    double phot_pack_start_Y; ///< the starting Y coordinate of photon package
    double phot_pack_start_Z; ///< the starting Z coordinate of photon package
};



class Photon_source
{
public:
    Photon_source();

	
	/** 
    *   @brief   Sets wave length of photon packet.   
    *  
    *   @param   double length - wave length of photon packet
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
    *   @brief   Sets Stokes I parameter for photon packet.   
    *  
    *   @param   double I_par - Stokes I parameter for photon packet
    *
	*   @return  void
    */  
    void set_stokes_I_param(double I_par);
    
	
	/** 
    *   @brief   Gets Stokes I parameter of photon packet.   
    *  
	*   @return  double - Stokes I parameter of photon packet
    */  
	double get_stokes_I_param();
	
	
	/** 
    *   @brief   Sets Stokes Q parameter for photon packet.   
    *  
    *   @param   double Q_par - Stokes Q parameter for photon packet
    *
	*   @return  void
    */  
    void set_stokes_Q_param(double Q_par);
    
	
	/** 
    *   @brief   Gets Stokes Q parameter of photon packet.   
    *  
	*   @return  double - Stokes Q parameter of photon packet
    */  
	double get_stokes_Q_param();
    
	
	/** 
    *   @brief   Sets Stokes U parameter for photon packet.   
    *  
    *   @param   double U_par - Stokes U parameter for photon packet
    *
	*   @return  void
    */  
	void set_stokes_U_param(double U_par);
    
	
	/** 
    *   @brief   Gets Stokes U parameter of photon packet.   
    *  
	*   @return  double - Stokes U parameter of photon packet
    */  
	double get_stokes_U_param();
    
	
	/** 
    *   @brief   Sets Stokes V parameter for photon packet.   
    *  
    *   @param   double V_par - Stokes V parameter for photon packet
    *
	*   @return  void
    */  
	void set_stokes_V_param(double V_par);
    
	
	/** 
    *   @brief   Gets Stokes V parameter of photon packet.   
    *  
	*   @return  double - Stokes V parameter of photon packet
    */  
	double get_stokes_V_param(); 
    
	
	/** 
    *   @brief   Sets photon weight for photon packet.   
    *  
    *   @param   double phot_weight - photon weight for photon packet.
    *
	*   @return  void
    */  
    void set_phot_weight(double phot_weight);

    
	/** 
    *   @brief   Generates polar and azimuth angles to choose surface
    *            point of photon packet emition.  
    *
	*   @return  void
    */  
    void generate_PHI_THET_surf();

    
	/** 
    *   @brief   Emmits photon packet.  
    *
	*   @return  emit_phot_param - structure that contains photon
	*            packet's parameters - Stokes parameters, wave length,
	*            azimuth and polar angles, photon weight, Cartesian
    *            coordinates (x,y,z) of photon package starting point. 
    */  
    emit_phot_param emit_phot();


    /** 
    *   @brief   Initialize photon source.  
    *
	*   @return  void
    */  
    void initialize_photon_source(star_coordinates *star_coord);


    /**
    *   @brief  Generates emition angles and emition point coordinates.
    * 
    *   @return void
    */
    void gener_emit_angles_point();


    /**
    *   @brief  Calculates point of photon packet emition.
    * 
    *   @return void
    */
    void calc_emition_point();


    /**
    *   @brief  Calculates unit vector of normal.
    * 
    *   @return void
    */
    void calc_normal_unit_vect();


    /**
    *   @brief  Generates angle delta_THET (between normal and emition vector)
    *           and alfa - angle of rotation about normal.
    * 
    *   @return void
    */
    void gener_delta_THET_alfa();


    /**
    *   @brief  Calculates vector of emition 0 (not final emition vector).
    * 
    *   @return void
    */
    void calc_emit_0_vect();


    /**
    *   @brief  Calculates vector of emition 1 (final emition vector).
    * 
    *   @return void
    */
    void calc_emit_1_vect();


    /**
    *   @brief  Calculates emition angles THET_1 and PHI_1.
    * 
    *   @return void
    */
    void calc_emition_angles();

    
private:

    double wave_length;      ///< wave length of photon packet
    
    double stokes_I_param;   ///< Stokes I parameter of photon packet
    double stokes_Q_param;   ///< Stokes Q parameter of photon packet
    double stokes_U_param;   ///< Stokes U parameter of photon packet
    double stokes_V_param;   ///< Stokes V parameter of photon packet
    
    double PHI_1;            ///< azimuth angle of photon packet
    double THET_1;           ///< polar angle of photon packet

    double photon_weight;    ///< quantity of photons in one photon packet

    double phot_pack_X;      ///< the starting X coordinate of photon package
    double phot_pack_Y;      ///< the starting Y coordinate of photon package
    double phot_pack_Z;      ///< the starting Z coordinate of photon package

    double star_center_X;    ///< X coordinate of star center
    double star_center_Y;    ///< Y coordinate of star center
    double star_center_Z;    ///< Z coordinate of star center

    double star_radius;      ///< radius of the star
    
    double PHI_surf;         ///< azimuth angle for choosing emition point on star surface
    double THET_surf;        ///< polar angle for choosing emition point on star surface

    double un_norm_vect_x;   ///< x parameter of unit vector of surface normal
    double un_norm_vect_y;   ///< y parameter of unit vector of surface normal
    double un_norm_vect_z;   ///< z parameter of unit vector of surface normal

    double un_em_0_vect_x;   ///< x parameter of emition vector 0 (that is not final emition vector)
    double un_em_0_vect_y;   ///< y parameter of emition vector 0 (that is not final emition vector)
    double un_em_0_vect_z;   ///< z parameter of emition vector 0 (that is not final emition vector)

    double un_em_1_vect_x;   ///< x parameter of emition vector 1 (that is final emition vector)
    double un_em_1_vect_y;   ///< y parameter of emition vector 1 (that is final emition vector)
    double un_em_1_vect_z;   ///< z parameter of emition vector 1 (that is final emition vector)

    double delta_THET;       ///< angle between normal vector and emition vector
    double alfa;             ///< angle of rotation of emition vector about normal vector
};


#endif