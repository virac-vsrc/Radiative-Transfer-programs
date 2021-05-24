/**
 * @file      angles_montecarlo.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Angles_montecarlo class functions.
 */

#include "angles_montecarlo.h"

using namespace std;

Angles_montecarlo::Angles_montecarlo()
{
    srand(time(NULL));
    init_Fxx_sphere();
    calc_decades();
}


void Angles_montecarlo::init_Fxx_sphere()
{
    ifstream infile1;
    
    infile1.open("F11F22F33F44F12F34");

    if(infile1.fail())
    {
        cerr << "opening File failed" << "\n";
    //    return 0;
    }

    int i = 0;

    while(!infile1.eof())
    {
        infile1 >> THET[i] >> F11[i] >> F22[i] >> F33[i] >> F44[i] >> F12[i] >> F34[i];
        i++;
    }

    infile1.close();
}


double Angles_montecarlo::func_2(double x)
{
    return x + F12[bigThetForFxx] / F11[bigThetForFxx] *
           (q_inc / 2.0 * sin(2*x) + u_inc/2.0*(cos(2*x)-1)) - 2.0*PI*eta; 
}


double Angles_montecarlo::bisection(double Epsilon) // limits from -0.01 to 6.3, result from 0 to 2*PI
{
    double a = 0;
    double b = 2*PI;
    double c = 0;
    
    int counter = 1;
    //eta = (rand() % 1000000) / 1000000.0;
    eta = (rand() % 1000) / 1000.0;

    while((b - a) >= Epsilon)
    {
        // Finding the middle point
        c = (a+b) / 2.0;

        // chose the limits of interval
        // the result of multiplication of limits should be negative
        if(func_2(c)*func_2(a) > 0)
            a = c;
        else
            b = c;

        counter++;
    }

    c = (a+b) / 2.0;
    
    //cout << " iteration number " << counter << "\n";
    psi = c;
    return c;
}


void Angles_montecarlo::set_iinc_qinc_uinc_vinc(double iinc, double qinc, double uinc, double vinc)
{
    i_inc = iinc;
    q_inc = qinc;
    u_inc = uinc;
    v_inc = vinc;
}


double Angles_montecarlo::get_q_sca()
{
    return q_sca;
}


double Angles_montecarlo::get_u_sca()
{
    return u_sca;
}


double Angles_montecarlo::get_v_sca()
{
    return v_sca;
}


void Angles_montecarlo::set_THET_1_PHI_1(double thet_1, double phi_1)
{
    THET_1 = thet_1;
    PHI_1 = phi_1;
}


double Angles_montecarlo::calc_THET_2()
{
    double val_for_acos = cos(THET_1)*cos(bigThet) - sin(THET_1)*sin(bigThet)*cos(psi);
    
    if(val_for_acos < -1)
    {
        val_for_acos = -1;
    }
    else if(val_for_acos > 1)
    {
        val_for_acos = 1;
    }

    THET_2 = acos(val_for_acos);
    return THET_2;
}


double Angles_montecarlo::calc_PHI_2() 
{
    PHI_2 = atan2( sin(bigThet)*sin(psi)/sin(THET_2),  (sin(THET_1)*cos(bigThet) + cos(THET_1)*sin(bigThet)*cos(psi))/sin(THET_2) );

    PHI_2 += PHI_1;

    if (PHI_2 < 0)
    {
        PHI_2 += 2*PI;
    }
    else if(PHI_2 > 2*PI)
    {
        PHI_2 -= 2*PI;
    }

    return PHI_2; 
}



double Angles_montecarlo::calc_cos_2psi()
{
    /* cout <<"cos(2*psi)=" << cos(2.0*psi) << "\n";
    cout <<"cos(2*psi)=" << 1.0 - 2.0*sin(THET_2)*sin(THET_2)
                            /(sin(bigThet)*sin(bigThet))
                            *sin(PHI_2-PHI_1)*sin(PHI_2-PHI_1) << "\n";
    */
    cos_2psi = cos(2*psi);

    return cos_2psi;
}


double Angles_montecarlo::calc_sin_2psi()
{
    /* cout << "sin(2*psi)=" << sin(2.0*psi) << "\n";
    cout << "sin(2*psi)=" << 2.0*sin(THET_2)*sin(PHI_2-PHI_1)
                             /(sin(bigThet)*sin(bigThet))
                             *(-sin(THET_1)*cos(THET_2) +
                             cos(THET_1)*sin(THET_2)*cos(PHI_2-PHI_1)) << "\n";
    */
    sin_2psi = sin(2*psi);

    return sin_2psi;
}


double Angles_montecarlo::calc_cos_2chi()
{
    cos_2chi = 1.0 - 2.0*sin(THET_1)*sin(THET_1)
                       *sin(PHI_2-PHI_1)*sin(PHI_2-PHI_1)
                       /(sin(bigThet)*sin(bigThet));

    return cos_2chi;
}

double Angles_montecarlo::calc_sin_2chi()
{
    sin_2chi = 2.0*sin(THET_1)*sin(PHI_2-PHI_1) / (sin(bigThet)*sin(bigThet))
               *(sin(THET_1)*cos(THET_2)*cos(PHI_2-PHI_1) - cos(THET_1)*sin(THET_2));
    
    return sin_2chi;
}


void Angles_montecarlo::calc_stokes_scat_par()
{
    calc_cos_2chi();
    calc_sin_2chi();
    calc_cos_2psi();
    calc_sin_2psi();

    
    stokes_vect[0] = i_inc;
    stokes_vect[1] = q_inc;
    stokes_vect[2] = u_inc;
    stokes_vect[3] = v_inc;


    matr_for_stokes_calc[0][0] = F11[bigThetForFxx];
    matr_for_stokes_calc[0][1] = F12[bigThetForFxx]*cos_2psi;
    matr_for_stokes_calc[0][2] = -F12[bigThetForFxx]*sin_2psi;
    matr_for_stokes_calc[0][3] = 0;

    matr_for_stokes_calc[1][0] = F12[bigThetForFxx]*cos_2chi;
    matr_for_stokes_calc[1][1] = F22[bigThetForFxx]*cos_2chi*cos_2psi - F33[bigThetForFxx]*sin_2chi*sin_2psi;
    matr_for_stokes_calc[1][2] = -F22[bigThetForFxx]*cos_2chi*sin_2psi - F33[bigThetForFxx]*sin_2chi*cos_2psi;
    matr_for_stokes_calc[1][3] = -F34[bigThetForFxx]*sin_2chi;
    
    matr_for_stokes_calc[2][0] = F12[bigThetForFxx]*sin_2chi;
    matr_for_stokes_calc[2][1] = F22[bigThetForFxx]*sin_2chi*cos_2psi + F33[bigThetForFxx]*cos_2chi*sin_2psi;
    matr_for_stokes_calc[2][2] = -F22[bigThetForFxx]*sin_2chi*sin_2psi + F33[bigThetForFxx]*cos_2chi*cos_2psi;
    matr_for_stokes_calc[2][3] = F34[bigThetForFxx]*cos_2chi;
    
    matr_for_stokes_calc[3][0] = 0;
    matr_for_stokes_calc[3][1] = -F34[bigThetForFxx]*sin_2psi;
    matr_for_stokes_calc[3][2] = -F34[bigThetForFxx]*cos_2psi;
    matr_for_stokes_calc[3][3] = F44[bigThetForFxx];


    double vect[4];
    vect[0] = 0;
    vect[1] = 0;
    vect[2] = 0;
    vect[3] = 0;

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            vect[i] += matr_for_stokes_calc[i][j]*stokes_vect[j];
        }
    }

    i_sca = vect[0] / vect[0];
    q_sca = vect[1] / vect[0];
    u_sca = vect[2] / vect[0];
    v_sca = vect[3] / vect[0];

    double sqrt_of_vect_sca = sqrt(q_sca*q_sca + u_sca*u_sca + v_sca*v_sca);

    if( 1 < sqrt_of_vect_sca)
    {
        q_sca = q_sca / sqrt_of_vect_sca;
        u_sca = u_sca / sqrt_of_vect_sca;
        v_sca = v_sca / sqrt_of_vect_sca;
    }

    sqrt_of_vect_sca = sqrt(q_sca*q_sca + u_sca*u_sca + v_sca*v_sca);

    if(i_sca + 0.00001 < sqrt_of_vect_sca )
    {
        cout << "i_sca = " << i_sca << "\n";
        cout << "q_sca = " << q_sca << "\n";
        cout << "u_sca = " << u_sca << "\n";
        cout << "v_sca = " << v_sca << "\n";

        cout << "bigThet = " << bigThet << "\n";
        cout << "THET_1 = " << THET_1 << "\n";
        cout << "PHI_1 = " << PHI_1 << "\n";
        cout << "THET_2 = " << THET_2 << "\n";
        cout << "PHI_2 = " << PHI_2 << "\n";
        cout << "psi = " << psi << "\n";
    }
}


double Angles_montecarlo::calc_bigThet()
{
    //xi = (rand() % 1000000) / 1000000.0;
    xi = (rand() % 1000) / 1000.0;
    //cout << "xi = " << xi << "\n";
    int i = 0;
    bigThet = 0;
    bigThetDeg = -1;
    double sum = 0;
    
    for(i = 0; i < 18; i++)
    {
        if(sumOfDecadesDeg[i] > xi * 2) // chose the starting integral
            break;
    }
    
    //cout << "we are using sumOfDecadesDeg[" << i <<"] =" << sumOfDecadesDeg[i] << "\n";

    if( i == 0) // BigThet is smaller than 10 degree
    {
        sum = 0;
        if((xi * 2) < 0.0001) // if xi is less than 0.25 then BigThet = 0 and out from cycle
        {
           // cout << "bigThet = " << bigThet << " out of cycle" << "\n";
        }
        else  // BigThet is bigger than zerro
        {
            for(int j = 0; j < 41; j++)
            {
                if(j == 0 || j == 40) // last number
                {
                    sum += PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
                }
                else if(j % 2 == 0) // even number
                {
                    sum += 2*PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
                }
                else // odd number
                {
                    sum += 4*PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);   
                }


                if(j % 2 == 0) // even number
                {
                    if((xi * 2) < sum) // xi is less than integral summ stop the cycle
                    {
                      //  cout << "xi is less than sum = " << sum << "\n";
                        sum -= PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
                        // the end for Simpson method
                        bigThet = j / 4.0 / 180.0 * PI;
                     //   cout << "sum = " << sum << "\n";
                        break; // out from "for" cycle
                    }
                }
            }
        }

    }
    else // i is bigger than 0
    {
        sum = sumOfDecadesDeg[i - 1] - PI/180.0/4.0/3.0 * F11[40*i]*sin(40*i*PI/180.0/4.0); 
        for(int j = 40 * i; j < 40*i + 41; j++)
        {
            if(j == (40*i + 40)) // last number
            {
                sum += PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
            }
            else if(j % 2 == 0) // even number
            {
                sum += 2*PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
            }
            else // odd number
            {
                sum += 4*PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);   
            }

            if(j % 2 == 0) // even number
            {
                if((xi*2) < sum)
                {
                  //  cout << "xi is less than sum = " << sum << "\n";
                    sum -= PI/180.0/4.0/3.0 * F11[j]*sin(j*PI/180.0/4.0);
                    // the end of Simpson method
                    bigThet = j / 4.0 / 180.0 * PI;
                   // cout << "sum = " << sum << "\n";
                    break; // out from "for" cycle
                }
            }
        }
    }

    if(bigThet < 1e-10)
    {
        bigThet = 1e-10;
    }
    else if(PI - bigThet < 1e-10)
    {
        bigThet = PI - 1e-10;
    }

    bigThetDeg = bigThet / PI * 180;
    bigThetForFxx = 4 * bigThetDeg + 0.5; // 0.5 is for rounding number

    if(bigThetForFxx > 720)
    {
        bigThetForFxx = 720;
    }

   // cout << "BigThetDeg is " << bigThetDeg << " deg. " << "\n";
   // cout << "bigThet is " << bigThet << " rad" << "\n";
    return bigThet;
}


void Angles_montecarlo::calc_next_decade(double previous_dec, int point_numb)
{
    sumOfDecadesDeg[point_numb / 40] = previous_dec + PI/180.0/4.0/3.0
                                        *F11[point_numb]*sin(point_numb*PI/180.0/4.0);

    for(int i = point_numb + 1; i < point_numb + 41; i++)
    {
        if(i == (point_numb + 40) ) // last number
        {
            sumOfDecadesDeg[point_numb / 40] += PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);
        }
        else if(i % 2 == 0) // even number
        {
            sumOfDecadesDeg[point_numb / 40] += 2*PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);
        }
        else // odd number
        {
            sumOfDecadesDeg[point_numb / 40] += 4*PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);
        }
    }
}


void Angles_montecarlo::calc_decades()
{
    sumOfDecadesDeg[0] = 0; // angle 10 degrees
    int i;
    for(i = 0; i < 41; i++)
    {
        if(i == 0 || i == 40) // first or last number
        {
            sumOfDecadesDeg[0] += PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);
        }
        else if(i % 2 == 0) // even number
        {
            sumOfDecadesDeg[0] += 2*PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);
        }
        else // odd number
        {
            sumOfDecadesDeg[0] += 4*PI/180.0/4.0/3.0 * F11[i]*sin(i*PI/180.0/4.0);   
        }
    }

    for(int j = 0; j < 18; j++)
    {
        calc_next_decade(sumOfDecadesDeg[j], (j+1)*40);
        //cout << (j+1)*10 << "deg sumOfDecadesDeg[" << j << "] = " << sumOfDecadesDeg[j] << "\n";
    }
}



float Angles_montecarlo::get_F11_val()
{
    return F11[Thet_peel_Fxx];
}



float Angles_montecarlo::get_F12_val()
{
    return F12[Thet_peel_Fxx];
}



void Angles_montecarlo::set_angles_for_peel(angles_for_peel angles_for_p)
{
    Thet_peel = angles_for_p.Thet_peeling;           
    Thet_tel = angles_for_p.Thet_teles;            
    Phi_tel = angles_for_p.Phi_teles;             
    Thet_inc = angles_for_p.Thet_incid;             
    Phi_inc = angles_for_p.Phi_incid;   

    Thet_peel_Fxx = (int)(Thet_peel / PI * 720 + 0.5); // +0.5 is for rounding number

    if(Thet_peel_Fxx > 720)
    {
        Thet_peel_Fxx = 720;
    }
}



double Angles_montecarlo::get_q_tel()
{
    return q_tel;
}



double Angles_montecarlo::get_u_tel()
{
    return u_tel;
}


   
double Angles_montecarlo::get_v_tel()
{
    return v_tel;
}



double Angles_montecarlo::calc_cos_2psi_tel()
{
    cos_2psi_tel = 1 - 2*sin(Thet_tel)*sin(Thet_tel)/(sin(Thet_peel)*sin(Thet_peel))
                   *sin(Phi_tel-Phi_inc)*sin(Phi_tel-Phi_inc);

    return cos_2psi_tel;
}
	


double Angles_montecarlo::calc_sin_2psi_tel()
{
    sin_2psi_tel = 2*sin(Thet_tel)*sin(Phi_tel-Phi_inc) / (sin(Thet_peel)*sin(Thet_peel))
                   *(cos(Thet_inc)*sin(Thet_tel)*cos(Phi_tel-Phi_inc) - sin(Thet_inc)*cos(Thet_tel));

    return sin_2psi_tel;
}

	
	
double Angles_montecarlo::calc_cos_2chi_tel()
{
    cos_2chi_tel = 1 - 2*sin(Thet_inc)*sin(Thet_inc)*sin(Phi_tel-Phi_inc)*sin(Phi_tel-Phi_inc)
                       / (sin(Thet_peel)*sin(Thet_peel));
    
    return cos_2chi_tel;
}
	
	
	
double Angles_montecarlo::calc_sin_2chi_tel()
{
    sin_2chi_tel = 2*sin(Thet_inc)*sin(Phi_tel-Phi_inc)/(sin(Thet_peel)*sin(Thet_peel))
                   *(sin(Thet_inc)*cos(Thet_tel)*cos(Phi_tel-Phi_inc)-cos(Thet_inc)*sin(Thet_tel));
    
    return sin_2chi_tel;
}


    
void Angles_montecarlo::calc_stokes_tel_par()
{
    calc_cos_2chi_tel();
    calc_sin_2chi_tel();
    calc_cos_2psi_tel();
    calc_sin_2psi_tel();


    stokes_vect_tel[0] = i_inc;
    stokes_vect_tel[1] = q_inc;
    stokes_vect_tel[2] = u_inc;
    stokes_vect_tel[3] = v_inc;


    matr_for_stokes_calc_tel[0][0] = F11[Thet_peel_Fxx];
    matr_for_stokes_calc_tel[0][1] = F12[Thet_peel_Fxx]*cos_2psi_tel;
    matr_for_stokes_calc_tel[0][2] = -F12[Thet_peel_Fxx]*sin_2psi_tel;
    matr_for_stokes_calc_tel[0][3] = 0;

    matr_for_stokes_calc_tel[1][0] = F12[Thet_peel_Fxx]*cos_2chi_tel;
    matr_for_stokes_calc_tel[1][1] = F22[Thet_peel_Fxx]*cos_2chi_tel*cos_2psi_tel - F33[Thet_peel_Fxx]*sin_2chi_tel*sin_2psi_tel;
    matr_for_stokes_calc_tel[1][2] = -F22[Thet_peel_Fxx]*cos_2chi_tel*sin_2psi_tel - F33[Thet_peel_Fxx]*sin_2chi_tel*cos_2psi_tel;
    matr_for_stokes_calc_tel[1][3] = -F34[Thet_peel_Fxx]*sin_2chi_tel;
    
    matr_for_stokes_calc_tel[2][0] = F12[Thet_peel_Fxx]*sin_2chi_tel;
    matr_for_stokes_calc_tel[2][1] = F22[Thet_peel_Fxx]*sin_2chi_tel*cos_2psi_tel + F33[Thet_peel_Fxx]*cos_2chi_tel*sin_2psi_tel;
    matr_for_stokes_calc_tel[2][2] = -F22[Thet_peel_Fxx]*sin_2chi_tel*sin_2psi_tel + F33[Thet_peel_Fxx]*cos_2chi_tel*cos_2psi_tel;
    matr_for_stokes_calc_tel[2][3] = F34[Thet_peel_Fxx]*cos_2chi_tel;
    
    matr_for_stokes_calc_tel[3][0] = 0;
    matr_for_stokes_calc_tel[3][1] = -F34[Thet_peel_Fxx]*sin_2psi_tel;
    matr_for_stokes_calc_tel[3][2] = -F34[Thet_peel_Fxx]*cos_2psi_tel;
    matr_for_stokes_calc_tel[3][3] = F44[Thet_peel_Fxx];


    double vect[4];
    vect[0] = 0;
    vect[1] = 0;
    vect[2] = 0;
    vect[3] = 0;

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            vect[i] += matr_for_stokes_calc_tel[i][j]*stokes_vect_tel[j];
        }
    }

    i_tel = vect[0] / vect[0];
    q_tel = vect[1] / vect[0];
    u_tel = vect[2] / vect[0];
    v_tel = vect[3] / vect[0];

    double sqrt_of_vect_tel = sqrt(q_tel*q_tel + u_tel*u_tel + v_tel*v_tel);

    if( 1 < sqrt_of_vect_tel)
    {
        q_tel = q_tel / sqrt_of_vect_tel;
        u_tel = u_tel / sqrt_of_vect_tel;
        v_tel = v_tel / sqrt_of_vect_tel;
    }

    sqrt_of_vect_tel = sqrt(q_tel*q_tel + u_tel*u_tel + v_tel*v_tel);

    if(i_tel + 0.00001 < sqrt_of_vect_tel)
    {
        cout << "i_tel = " << i_tel << "\n";
        cout << "q_tel = " << q_tel << "\n";
        cout << "u_tel = " << u_tel << "\n";
        cout << "v_tel = " << v_tel << "\n";

        cout << "Thet_peel = " << Thet_peel << "\n";
        cout << "Thet_inc = " << Thet_inc << "\n";
        cout << "Phi_inc = " << Phi_inc << "\n";
        cout << "Thet_tel = " << Thet_tel << "\n";
        cout << "Phi_tel = " << Phi_tel << "\n";
        cout << "bigThet = " << bigThet << "\n";
    }
}