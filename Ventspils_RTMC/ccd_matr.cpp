/**
 * @file      ccd_matr.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of CCD_matr class functions.
 */

#include "ccd_matr.h"


CCD_matr::CCD_matr()
{
    for(int i = 0; i < ARRAY_X*CCD_PIXELS; i++)
    {
        for(int j = 0; j < ARRAY_Y*CCD_PIXELS; j++)
        {
            matrix[i][j].photon_weight = 0;
            matrix[i][j].I_stokes = 0;
            matrix[i][j].q_stokes = 0;
            matrix[i][j].u_stokes = 0;
            matrix[i][j].v_stokes = 0;
        }
    }
}


void CCD_matr::add_cell_phot_weight(int x_cor, int y_cor,long double phot_weight)
{
    lock_guard<mutex> locker_1(mu_scattering);

    matrix[x_cor][y_cor].photon_weight += phot_weight;
}


void CCD_matr::add_cell_I_sto(int x_cor, int y_cor,long double I_par)
{
    lock_guard<mutex> locker_8(mu_i);
    matrix[x_cor][y_cor].I_stokes += I_par;
}


void CCD_matr::add_cell_q_sto(int x_cor, int y_cor,long double q_par)
{
    lock_guard<mutex> locker_9(mu_q);
    matrix[x_cor][y_cor].q_stokes += q_par;
}


void CCD_matr::add_cell_u_sto(int x_cor, int y_cor,long double u_par)
{
    lock_guard<mutex> locker_10(mu_u);
    matrix[x_cor][y_cor].u_stokes += u_par;
}


void CCD_matr::add_cell_v_sto(int x_cor, int y_cor,long double v_par)
{
    lock_guard<mutex> locker_11(mu_v);
    matrix[x_cor][y_cor].v_stokes += v_par;
}


long double CCD_matr::get_cell_phot_weight(int x_cor, int y_cor)
{
    return matrix[x_cor][y_cor].photon_weight;
}


long double CCD_matr::get_cell_I_sto(int x_cor, int y_cor)
{
    return matrix[x_cor][y_cor].I_stokes;
}


long double CCD_matr::get_cell_q_sto(int x_cor, int y_cor)
{
    return matrix[x_cor][y_cor].q_stokes;
}


long double CCD_matr::get_cell_u_sto(int x_cor, int y_cor)
{
    return matrix[x_cor][y_cor].u_stokes;
}


long double CCD_matr::get_cell_v_sto(int x_cor, int y_cor)
{
    return matrix[x_cor][y_cor].v_stokes;
}




void CCD_matr::generate_gnu_txt_file()
{
    ofstream data_gnu_I;
    data_gnu_I.open("data_gnu_I.txt");

    for(int i = 0; i < ARRAY_Y*CCD_PIXELS; i++)
    {
        for(int j = 0; j < ARRAY_X*CCD_PIXELS; j++)
        {
            data_gnu_I << matrix[j][i].photon_weight << " "; // CCD matrix cells for Stokes vector I parameter
        }
        data_gnu_I << "\n";
    }
    data_gnu_I.close();


    ofstream data_gnu_Q;
    data_gnu_Q.open("data_gnu_Q.txt");

    for(int i = 0; i < ARRAY_Y*CCD_PIXELS; i++)
    {
        for(int j = 0; j < ARRAY_X*CCD_PIXELS; j++)
        {
            data_gnu_Q << matrix[j][i].q_stokes << " "; // CCD matrix cells for Stokes vector Q parameter
        }
        data_gnu_Q << "\n";
    }
    data_gnu_Q.close();


    ofstream data_gnu_U;
    data_gnu_U.open("data_gnu_U.txt");

    for(int i = 0; i < ARRAY_Y*CCD_PIXELS; i++)
    {
        for(int j = 0; j < ARRAY_X*CCD_PIXELS; j++)
        {
            data_gnu_U << matrix[j][i].u_stokes << " "; // CCD matrix cells for Stokes vector U parameter
        }
        data_gnu_U << "\n";
    }
    data_gnu_U.close();


    ofstream data_gnu_V;
    data_gnu_V.open("data_gnu_V.txt");

    for(int i = 0; i < ARRAY_Y*CCD_PIXELS; i++)
    {
        for(int j = 0; j < ARRAY_X*CCD_PIXELS; j++)
        {
            data_gnu_V << matrix[j][i].v_stokes << " "; // CCD matrix cells for Stokes vector V parameter
        }
        data_gnu_V << "\n";
    }
    data_gnu_V.close();
}