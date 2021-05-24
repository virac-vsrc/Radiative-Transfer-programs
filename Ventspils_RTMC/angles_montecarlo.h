/**
 * @file      angles_montecarlo.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class to calculate the polar, azimuthal and other angles 
 *	          for photon packet.
 *
 * @section   DESCRIPTION
 *
 *            This class is essential to calculate photon packet travelling
 *			  angles after each scattering. The angles are calculating using
 *			  pseido random number generator. To calculate scattering angle
 *			  bigThet the integral equation need to be solved, comparing the
 *			  integral result with pseido random number xi. To calculate 
 *            angle psi the transedental equation is calculated using 
 *            bisection method, the pseido random number eta is used in 
 *            equation. After calculating angle bigThet and psi the polar
 *            scattering angle THET_2 could be calculated using tigonometric
 *            equationg. After calculating angle THET_2 the azimuth angle
 *            PHI_2 could be calculated. Class is used to calculate Stokes
 *            vectors for photon packet and for telescope.
 */
 
#ifndef ANGLES_MONTECARLO_H_INCLUDED
#define ANGLES_MONTECARLO_H_INCLUDED

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include "dust_cloud.h"
#include "photon_packet.h"


class Angles_montecarlo
{
    public:
    Angles_montecarlo();

	
	/** 
    *   @brief  Initializes arrays with angle and Stokes scattering matrix elements   
    *  
    *   @return  void
    */  
    void init_Fxx_sphere();

	
	/** 
    *   @brief  Calculates scattering matrix elements integrals for 10, 20 ... 180 degrees  
    *  
    *   @return  void
    */  
    void calc_decades();
	
	
	/** 
    *   @brief   Calculates decade for scattering matrix elements integral.   
    *  
    *   @param   double previous_dec - the result of integral of previous decade
	*   @param   int point_numb - the number of scattering matrix array element
    *
	*   @return  void
    */  
    void calc_next_decade(double previous_dec, int point_numb);
    
	
    /** 
    *   @brief   Solves the transedental equation to calculate angle psi.
    *            Angle psi is essential to calculate polar angle THET_2 
    *            and azimuth angle PHI_2. Psi is needed to rotate 
    *            Stokes vector coordinate system.
    *  
    *   @param   double Epsilon - the accuracy of equation calculation
    *
	*   @return  double - the value of psi
    */  
    double bisection(double Epsilon);

	
	/** 
    *   @brief   The equation for angle psi.   
    *  
    *   @param   double x - the equation root
    *
	*   @return  double - the result of equation (zerro if root is correct)
    */  
    double func_2(double x); // the second equation

	
	/** 
    *   @brief   Calculates the scattering angle bigThet.   
    *  
    *   @return  double - the bigThet angle in radians
    */  
    double calc_bigThet();

	
	/** 
    *   @brief   Sets I, Q, U and V Stokes parameters divided by I.  
    *   
    *   @param   double iinc - Stokes parameter I divided by I  = 1
    *   @param   double qinc - Stokes parameter Q divided by I
	*   @param   double uinc - Stokes parameter U divided by I
    *   @param   double vinc - Stokes parameter V divided by I
    *
	*   @return  void
    */  
    void set_iinc_qinc_uinc_vinc(double iinc, double qinc, double uinc, double vinc);
    

    /** 
    *   @brief   Gets calculated Stokes parameter Q divided by I, for scattered photon packet.  
    *
	*   @return  double - calculated Stokes parameter Q divided by I
    */  
    double get_q_sca();


    /** 
    *   @brief   Gets calculated Stokes parameter U divided by I, for scattered photon packet.  
    *
	*   @return  double - calculated Stokes parameter U divided by I
    */  
    double get_u_sca();


    /** 
    *   @brief   Gets calculated Stokes parameter V divided by I, for scattered photon packet.  
    *
	*   @return  double - calculated Stokes parameter V divided by I
    */  
    double get_v_sca();
	

    /** 
    *   @brief   Sets incident polar angle THET_1 and incident azimuth angle PHI_1 of photon packet.
    *  
    *   @param   double thet_1 - incident polar angle
	*   @param   double phi_1 - incident azimuth angle
    *
	*   @return  void
    */  
    void set_THET_1_PHI_1(double thet_1, double phi_1);
    

	/** 
    *   @brief   Calculates scattering polar angle THET_2 of photon packet.
    *  
    *   @return  double - polar angle THET_2
    */  
    double calc_THET_2();

	
	/** 
    *   @brief   Calculates scattering azimuth angle PHI_2 of photon packet.
    *  
    *   @return  double - azimuth angle PHI_2
    */  
    double calc_PHI_2();

	
	/** 
    *   @brief   Calculates cosine of 2*psi with two different methods. 
    *  
    *   @return  double - cosine of 2*psi
    */  
    double calc_cos_2psi();
	
	
	/** 
    *   @brief   Calculates sine of 2*psi with two different methods.
    *  
    *   @return  double - sine of 2*psi
    */  
    double calc_sin_2psi();

	
	/** 
    *   @brief   Calculates cosine of 2*chi. 
    *  
    *   @return  double - cosine of 2*chi
    */  
    double calc_cos_2chi();
	
	
	/** 
    *   @brief   Calculates sine of 2*chi
    *  
    *   @return  double - sine of 2*chi
    */  
    double calc_sin_2chi();


    /** 
    *   @brief   Calculates Stokes parameters Q, U ,V divided by I, for scattered photon packet.
    *  
    *   @return  void
    */  
    void calc_stokes_scat_par();


    /** 
    *   @brief   Gets F11 value for Thet_val angle.
    *  
    *   @return  float - F11 value.
    */ 
    float get_F11_val();


    /** 
    *   @brief   Gets F12 value for Thet_val angle.
    *  
    *   @return  float - F12 value.
    */ 
    float get_F12_val();


    /** 
    *   @brief   Sets angles for peeling off method.
    *  
    *   @param   angles_for_peel angles_for_p - structure of angles for peeling off method:
    *                                           Scattering angle in telescope direction;
    *                                           Polar angle in telescope direction
    *                                           Azimuth angle in telescope direction;
    *                                           Polar angle of incident photon package;
    *                                           Azimuth angle of incident photon package.
    * 
    *   @return  void
    */ 
    void set_angles_for_peel(angles_for_peel angles_for_p);


    /** 
    *   @brief   Gets calculated Stokes parameter Q divided by I, for telescope.  
    *
	*   @return  double - calculated Stokes parameter Q divided by I
    */  
    double get_q_tel();


    /** 
    *   @brief   Gets calculated Stokes parameter U divided by I, for telescope.  
    *
	*   @return  double - calculated Stokes parameter U divided by I
    */  
    double get_u_tel();


    /** 
    *   @brief   Gets calculated Stokes parameter V divided by I, for telescope.  
    *
	*   @return  double - calculated Stokes parameter V divided by I
    */  
    double get_v_tel();


    /** 
    *   @brief   Calculates cosine of 2*psi_tel. 
    *  
    *   @return  double - cosine of 2*psi_tel
    */  
    double calc_cos_2psi_tel();
	
	
	/** 
    *   @brief   Calculates sine of 2*psi_tel.
    *  
    *   @return  double - sine of 2*psi_tel
    */  
    double calc_sin_2psi_tel();

	
	/** 
    *   @brief   Calculates cosine of 2*chi_tel
    *  
    *   @return  double - cosine of 2*chi_tel
    */  
    double calc_cos_2chi_tel();
	
	
	/** 
    *   @brief   Calculates sine of 2*chi_tel
    *  
    *   @return  double - sine of 2*chi_tel
    */  
    double calc_sin_2chi_tel();


    /** 
    *   @brief   Calculates Stokes parameters Q, U ,V divided by I, for telescope.
    *  
    *   @return  void
    */  
    void calc_stokes_tel_par();


    private:

    double xi;                  ///< random number from [0;1]
    double eta;                 ///< random number from [0;1]
    double psi;                 ///< angle that essential to calculate THET_2 and PHI_2, 
                                ///< and needed to rotate Stokes vector coordinate system
    double chi;                 ///< angle that is needed to rotate Stokes vector coordinate system
    double cos_2chi;            ///< calculated cosine of 2*chi
    double sin_2chi;            ///< calculated sine of 2*chi
    double cos_2psi;            ///< calculated cosine of 2*psi
    double sin_2psi;            ///< calculated sine of 2*psi
    double bigThet;             ///< scattering angle of photon packet rad
    double bigThetDeg;          ///< bigThet in degrees
    int bigThetForFxx;          ///< bigThet for F11 and F12 arrays

    double i_inc;               ///< Stokes parameter I_inc/I_inc = 1 of incident photon packet
    double q_inc;               ///< Stokes parameter Q_inc/I_inc of incident photon packet
    double u_inc;               ///< Stokes parameter U_inc/I_inc of incident photon packet
    double v_inc;               ///< Stokes parameter V_inc/I_inc of incident photon packet

    double i_sca;               ///< Stokes parameter I_inc/I_inc = 1 of scattered photon packet
    double q_sca;               ///< Stokes parameter Q_inc/I_inc of scattered photon packet
    double u_sca;               ///< Stokes parameter U_inc/I_inc of scattered photon packet
    double v_sca;               ///< Stokes parameter V_inc/I_inc of scattered photon packet
    
    double matr_for_stokes_calc[4][4];  ///< Matrix for Stokes parameters calculation
    double stokes_vect[4];              ///< Vector with Stokes parameters divided by I

    float THET[721];            ///< angle bigThetDeg from 0 to 180 degrees with step 0.25 degree
    float F11[721];             ///< F11 parameter for spheroidal particles with step 0.25 degree
    float F22[721];             ///< F22 parameter for spheroidal particles with step 0.25 degree
    float F33[721];             ///< F33 parameter for spheroidal particles with step 0.25 degree
    float F44[721];             ///< F44 parameter for spheroidal particles with step 0.25 degree
    float F12[721];             ///< F12 parameter for spheroidal particles with step 0.25 degree
    float F34[721];             ///< F34 parameter for spheroidal particles with step 0.25 degree

    double sumOfDecadesDeg[18]; ///< sum of integrals for 10, 20, 30 ... 180 degree
    double sum;                 ///< esential to calculate integrals for 10, 20, 30 ... 180 degree

    double THET_1;              ///< incident polar angle of photon packet [0;PI]  in radians
    double PHI_1;               ///< incident angle of incident photon packet [0;2*PI] in radians
    double THET_2;              ///< scattering polar angle of photon packet [0;PI] in radians
    double PHI_2;               ///< scattering azimuth angle of photon packet[0;2*PI] in radians


    //---for peeling off method ---//
    double Thet_peel;           ///< Scattering angle in telescope direction
    float Thet_tel;             ///< Polar angle in telescope direction
    float Phi_tel;              ///< Azimuth angle in telescope direction
    float Thet_inc;             ///< Polar angle of incident photon package
    float Phi_inc;              ///< Azimuth angle of incident photon package

    int Thet_peel_Fxx;          ///< Scattering angle in telescope direction for Fxx

    double i_tel;               ///< Stokes parameter I_tel / I_tel = 1 , for telescope
    double q_tel;               ///< normalized Stokes vector Q parameter, for telescope
    double u_tel;               ///< normalized Stokes vector U parameter, for telescope
    double v_tel;               ///< nomralized Stokes vector V parameter, for telescope

    double cos_2chi_tel;            ///< calculated cosine of 2*chi_tel
    double sin_2chi_tel;            ///< calculated sine of 2*chi_tel
    double cos_2psi_tel;            ///< calculated cosine of 2*psi_tel
    double sin_2psi_tel;            ///< calculated sine of 2*psi_tel

    double matr_for_stokes_calc_tel[4][4];  ///< Matrix for Stokes parameters calculation for telescpe
    double stokes_vect_tel[4];              ///< Vector with Stokes parameters divided by I for telescope
    //----------------------------//

};

#endif