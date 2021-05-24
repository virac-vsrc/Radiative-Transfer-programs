/**
 * @file      dust_cloud.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 * @date 
 * @brief     Description of Dust_cloud class functions.
 */


#include "dust_cloud.h"

using namespace std;


Dust_cloud::Dust_cloud(short cloud_Type)
{
    cloud_type = cloud_Type;
    switch(cloud_type)
    {
        case 1:
            cout << "You have chosen the sphere\n";          
            break;
        case 2:
            cout << "You have chosen assymetrical_1 shape \n";
        default: 
            break;

    }

    for(int i = 0; i < ARRAY_X; i++)
    {
        for(int j = 0; j < ARRAY_Y; j++)
        {
            for(int k = 0; k < ARRAY_Z; k++)
            {
                three_dim_array[i][j][k] = 0;
            }
        }
    }
}


Dust_cloud::~Dust_cloud()
{}


double Dust_cloud::sphere_X_func(double phi, double thet)
{
    return outer_rad*cos(phi)*sin(thet) + ARRAY_X/2;
}


double Dust_cloud::sphere_Y_func(double phi, double thet)
{
    return outer_rad*sin(phi)*sin(thet) + ARRAY_Y/2;
}


double Dust_cloud::sphere_Z_func(double thet)
{
    return outer_rad*cos(thet) + ARRAY_Z/2;
}


double Dust_cloud::assymetrical_1_X_func(double thet)
{
    return outer_rad*cos(thet) + ARRAY_X/2;
}
	
	
double Dust_cloud::assymetrical_1_Y_func(double phi, double thet)
{
    return outer_rad*(sin(phi)*sin(thet)-sin(thet/3)) + ARRAY_Y*3/4;
}
	
	
double Dust_cloud::assymetrical_1_Z_func(double phi, double thet)
{
    return outer_rad*(cos(phi)*sin(thet)+cos(thet/2) ) + ARRAY_Z/4;
}


void Dust_cloud::set_min_max_phi_thet(double min_Phi, double min_Thet, double max_Phi, double max_Thet)
{
    min_phi = min_Phi;
    min_thet = min_Thet;
    max_phi = max_Phi;
    max_thet = max_Thet;
}

void Dust_cloud::set_outer_rad(double out_Rad)
{
    outer_rad = out_Rad;
}


short Dust_cloud::get_array_value(short i, short j, short k)
{
    if(i >= 0  &&  i < ARRAY_X  &&  j >= 0  &&  j < ARRAY_Y  &&  k >= 0  &&  k < ARRAY_Z)
    {
        return three_dim_array[i][j][k];
    }    
    else
    {
        cout << "Wrong coordinates of an array.\n";
        return 9999;
    }
}


void Dust_cloud::make_cloud_surface()
{
    float step_phi = (max_phi - min_phi) / 1440.0; // azimuth angle
    float step_thet = (max_thet - min_thet) / 720.0; // polar angle

    double x = 0;
    double y = 0;
    double z = 0;

    switch(cloud_type)
    {
        case 1:
            cout << "The sphere surface is created\n";
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = sphere_X_func(i, j);
                    y = sphere_Y_func(i, j);
                    z = sphere_Z_func(j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }

            outer_rad++;
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = sphere_X_func(i, j);
                    y = sphere_Y_func(i, j);
                    z = sphere_Z_func(j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }

            outer_rad -= 2;
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = sphere_X_func(i, j);
                    y = sphere_Y_func(i, j);
                    z = sphere_Z_func(j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }
            outer_rad++;


            break;
        case 2:
            cout << "The assymetrical_1 surface is created\n";
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = assymetrical_1_X_func(j);
                    y = assymetrical_1_Y_func(i, j);
                    z = assymetrical_1_Z_func(i,j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }

            outer_rad++;
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = assymetrical_1_X_func(j);
                    y = assymetrical_1_Y_func(i, j);
                    z = assymetrical_1_Z_func(i, j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }

            outer_rad -= 2;
            for(float i = min_phi; i <= max_phi; i+= step_phi)
            {
                for(float j = min_thet; j <= max_thet; j+=step_thet)
                {
                    x = assymetrical_1_X_func(j);
                    y = assymetrical_1_Y_func(i, j);
                    z = assymetrical_1_Z_func(i, j);
                    three_dim_array[(int)x][(int)y][(int)z] = 1;
                }            
            }
            outer_rad++;


            break;
        default:
            break;
    }
}



void Dust_cloud::cut_off_sphere(star_coordinates *st_coord[star_quantity])
{
    short star_numb = 0;

    while(st_coord[star_numb] != NULL)
    {
        short x = 0;
        short y = 0;
        short z = 0;
    
        float phi = 2*PI; // azimuth angle for sphere 
        float thet = PI; // polar angle for sphere
        short inner_Radius_1 = (short) st_coord[star_numb]->cloud_radius; // inside of cloud is photon source and gas without dusts
        short inner_Radius_2 = (short) st_coord[star_numb]->cloud_radius + 1;
        //short inner_Radius_3 = (short) st_coord[star_numb]->cloud_radius - 1;

        for(float i = 0; i <= phi; i += phi / 720.0)
        {
            for(float j = 0; j <= thet; j += thet / 360.0)
            {
                x = (int)(inner_Radius_1 * cos(i) * sin(j) - 0.5) + st_coord[star_numb]->x_pos;
                y = (int)(inner_Radius_1 * sin(i) * sin(j) - 0.5) + st_coord[star_numb]->y_pos;
                z = (int)(inner_Radius_1 * cos(j) - 0.5) + st_coord[star_numb]->z_pos;
                three_dim_array[x][y][z] = 3 + star_numb;  
                x = (int)(inner_Radius_2 * cos(i) * sin(j) - 0.5) + st_coord[star_numb]->x_pos;
                y = (int)(inner_Radius_2 * sin(i) * sin(j) - 0.5) + st_coord[star_numb]->y_pos;
                z = (int)(inner_Radius_2 * cos(j) - 0.5) + st_coord[star_numb]->z_pos;
                three_dim_array[x][y][z] = 3 + star_numb;    
            }
        }

        // fill the Volume that cut off from dust cloud
        for(int i = 0; i < ARRAY_X; i++)
        {
            for(int j = 0; j < ARRAY_Y; j++)
            {
                for(int k = 0; k < ARRAY_Z; k++)
                {
                    // unset variables
                    cloud_UP = 0; 
                    cloud_DOWN = 0;
                    cloud_FRONT = 0;
                    cloud_BACK = 0;
                    cloud_RIGHT = 0;
                    cloud_LEFT = 0;
                    // fill an array elements that cut off from dust cloud
                    if(three_dim_array[i][j][k] == 0)
                    {
                        for(x = i+1; x < ARRAY_X; x++)
                        {
                            if( three_dim_array[x][j][k] == 3 + star_numb)
                            {
                                cloud_RIGHT = 1;
                            }
                        }
                    
                        for(x = i-1; x >= 0; x--)
                        {
                            if( three_dim_array[x][j][k] == 3 + star_numb)
                            {
                                cloud_LEFT = 1;
                            }
                        }

                        for(y = j+1; y < ARRAY_Y; y++)
                        {
                            if( three_dim_array[i][y][k] == 3 + star_numb)
                            {
                                cloud_FRONT = 1;
                            }
                        }

                        for(y = j-1; y >= 0; y--)
                        {
                            if( three_dim_array[i][y][k] == 3 + star_numb)
                            {
                                cloud_BACK = 1;
                            }
                        }

                        for(z = k+1; z < ARRAY_Z; z++)
                        {
                            if( three_dim_array[i][j][z] == 3 + star_numb)
                            {
                                cloud_UP = 1;
                            }
                        }

                        for(z = k-1; z >= 0; z--)
                        {
                            if( three_dim_array[i][j][z] == 3 + star_numb)
                            {
                                cloud_DOWN = 1;
                            }
                        }
                    }

                    if(cloud_UP == 1 && cloud_DOWN == 1 && cloud_FRONT == 1 && cloud_BACK == 1
                         && cloud_RIGHT == 1 && cloud_LEFT == 1 && three_dim_array[i][j][k]!= 3)
                    {
                        three_dim_array[i][j][k] = 40 + star_numb;
                    }
                }
            }
        }

        star_numb++;
    }
}



void Dust_cloud::fill_cloud()
{
// fill the cloud of dust Volume with dusts

    short x = 0;
    short y = 0;
    short z = 0;

    
    for(int i = 0; i < ARRAY_X; i++)
    {
        for(int j = 0; j < ARRAY_Y; j++)
        {
            for(int k = 0; k < ARRAY_Z; k++)
            {
                // unset variables
                cloud_UP = 0; 
                cloud_DOWN = 0;
                cloud_FRONT = 0;
                cloud_BACK = 0;
                cloud_RIGHT = 0;
                cloud_LEFT = 0;
                // fill an array elements that are inside the surface
                if(three_dim_array[i][j][k] == 0)
                {
                    for(x = i+1; x < ARRAY_X; x++)
                    {
                        if( three_dim_array[x][j][k] == 1)
                        {
                            cloud_RIGHT = 1;
                        }
                    }

                    for(x = i-1; x >= 0; x--)
                    {
                        if( three_dim_array[x][j][k] == 1)
                        {
                            cloud_LEFT = 1;
                        }
                    }

                    for(y = j+1; y < ARRAY_Y; y++)
                    {
                        if( three_dim_array[i][y][k] == 1)
                        {
                            cloud_FRONT = 1;
                        }
                    }

                    for(y = j-1; y >= 0; y--)
                    {
                        if( three_dim_array[i][y][k] == 1)
                        {
                            cloud_BACK = 1;
                        }
                    }

                    for(z = k+1; z < ARRAY_Z; z++)
                    {
                        if( three_dim_array[i][j][z] == 1)
                        {
                            cloud_UP = 1;
                        }
                    }

                    for(z = k-1; z >= 0; z--)
                    {
                        if( three_dim_array[i][j][z] == 1)
                        {
                            cloud_DOWN = 1;
                        }
                    }
                }

                if(cloud_UP == 1 && cloud_DOWN == 1 && cloud_FRONT == 1 && cloud_BACK == 1
                     && cloud_RIGHT == 1 && cloud_LEFT == 1 && three_dim_array[i][j][k]!= 1)
                {
                    three_dim_array[i][j][k] = 2; // to see which points are changed
                }
            }
        }
    }
    
    /*//////////////////////
    z = ARRAY_Z / 2;
    //for(z = ARRAY_Z/2; z > 0; z--)
    //{
        for( x = 0; x < ARRAY_X; x++)
        {
            for( y = 0; y < ARRAY_Y; y++)
            {
                //if(three_dim_array[x][y][z] == 3)
                //{
                //    cout << "4"; 
                //}
                //else
                //{
                    cout << three_dim_array[x][y][z]; // show the layers of 3 dimensional array that describes dust cloud
                //}
            
            }
            cout << "\n";
        } 
        //cout << "z = " << z << "\n";
    //}
    /*/////////////////// 
}