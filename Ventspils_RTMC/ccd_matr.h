/**
 * @file      ccd_matr.h
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     The class to simulate virtual CCD matrix.
 *
 * @section   DESCRIPTION
 *
 *            The class essential to collect all scattering data of photon
 *            packets. Virtual CCD matrix contains 2-dimensional array
 *            of structure that contains Stokes vector parameters and 
 *            photon weight parameter. After each photon packet 
 *            scattering the parameters are added to CCD matrix cell.
 *            After emitting and scattering photon packets the data of
 *            light scattering is collect. This data is saved in file. 
 *            Using Python language the data file is converted to image.
 *
 */

#ifndef CCD_MATR_H_INCLUDED
#define CCD_MATR_H_INCLUDED

#include <iostream>
#include <fstream>
#include <math.h>
#include "dust_cloud.h"
#include <mutex>

//#define CCD_PIXELS 1.25
#define CCD_PIXELS 2.5
//#define CCD_PIXELS 9               ///< how many pixels in one 3-dimensional array element. 9x9 = 81 pixel
//#define CCD_PIXELS 11            ///< how many pixels in one 3-dimensional array element. 11x11 = 121 pixel
//#define CCD_PIXELS 12            ///< how many pixels in one 3-dimensional array element. 12x12 = 144 pixel
//#define CCD_PIXELS 12.5          ///< how many pixels in one 3-dimensional array element. 12.5x12.5 = 156.25 pixel 
//#define CCD_PIXELS 25            ///< how many pixels in one 3-dimensional array element. 25x25 = 625 pixel
//#define CCD_PIXELS 50            ///< how many pixels in one 3-dimensional array element. 50x50 = 250 pixel

using namespace std;

struct one_cell                   ///< one cell of virtual CCD matrix
{
    long double photon_weight;    ///< photon weight of one cell
    long double I_stokes;         ///< I Stokes parameter of one cell
    long double q_stokes;         ///< Q Stokes parameter divided by I (for one cell)
    long double u_stokes;         ///< U Stokes parameter divided by I (for one cell)
    long double v_stokes;         ///< V Stokes parameter divided by I (for one cell)
};

class CCD_matr
{
    public:

    CCD_matr();

	
	/** 
    *   @brief   Adds photon weight parameter for one CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell
	*   @param   double phot_weight - the photon weight parameter
    *
	*   @return  void
    */  
    void add_cell_phot_weight(int x_cor, int y_cor,long double phot_weight);
    
	
	/** 
    *   @brief   Adds I Stokes parameter for one CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell  
	*   @param   double I_par - Stokes I parameter
    *
	*   @return  void
    */  
	void add_cell_I_sto(int x_cor, int y_cor,long double I_par);
    
	
	/** 
    *   @brief   Adds Q Stokes parameter divided by I parameter for one
	*            CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell 
	*   @param   double q_par - Q Stokes parameter divided by I 
	*            parameter
    *
	*   @return  void
    */  
	void add_cell_q_sto(int x_cor, int y_cor,long double q_par);
    
	
	/** 
    *   @brief   Adds U Stokes parameter divided by I parameter for one
	*            CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell  
	*   @param   double u_par - U Stokes parameter divided by I
	*            parameter
    *
	*   @return  void
    */  
	void add_cell_u_sto(int x_cor, int y_cor,long double u_par);
    
	
	/** 
    *   @brief   Adds V Stokes parameter divided by I parameter for one
    *            CCD matrix cell.	
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell 
	*   @param   double v_par - V Stokes parameter divided by I_par
	*            parameter
    *
	*   @return  void
    */  
	void add_cell_v_sto(int x_cor, int y_cor,long double v_par);

	
	/** 
    *   @brief   Gets photon weight of one CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell 	
    *
	*   @return  double - photon weight of one CCD matrix cell
    */  
    long double get_cell_phot_weight(int x_cor, int y_cor);
	
	
	/** 
    *   @brief   Gets Stokes I parameter of one CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell 	
    *
	*   @return  double - Stokes I parameter of one CCD matrix cell
    */  
    long double get_cell_I_sto(int x_cor, int y_cor);
    
	
	/** 
    *   @brief   Gets Stokes Q parameter divided by I parameter of one
	*            CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell	
    *
	*   @return  double - Stokes Q parameter divided by I parameter of
	*            one CCD matrix cell
    */  
	long double get_cell_q_sto(int x_cor, int y_cor);
    
	
	/** 
    *   @brief   Gets Stokes U parameter divided by I parameter of one
	*            CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell 
    *
	*   @return  double - Stokes U parameter divided by I parameter of
	*                     one CCD matrix cell.
    */  
	long double get_cell_u_sto(int x_cor, int y_cor);
    
	
	/** 
    *   @brief   Gets Stokes V parameter divided by I parameter of one
	*            CCD matrix cell.
    *  
    *   @param   short x_cor - x coordinate of CCD matrix cell
	*   @param   short y_cor - y coordinate of CCD matrix cell   	
    *
	*   @return  double - Stokes V parameter divided by I parameter of
	*            one CCD matrix cell.
    */  
	long double get_cell_v_sto(int x_cor, int y_cor);

	
	/** 
    *   @brief   Generates file with virtual CCD matrix data for GNUplot.
	*            It is possible to convert data file into image using 
	*            GNUplot scripts.
    *  
	*   @return  void
    */  
    void generate_gnu_txt_file();

    private:
    one_cell matrix[ (int) (ARRAY_X * CCD_PIXELS)][ (int) (ARRAY_Y * CCD_PIXELS)];  ///< 2-dimensional array - Virtual CCD matrix

    mutex mu_scattering;
    mutex mu_i;
    mutex mu_q;
    mutex mu_u;
    mutex mu_v;
};

#endif