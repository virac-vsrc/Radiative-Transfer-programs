/**
 * @file      main.cpp
 * @authors   Juris Freimanis and Romans Pezenkovs
 */

#include <iostream>
#include <iomanip>
#include "dust_cloud.h"
#include "sub_Volume.h"
#include "photon_source.h"
#include "photon_packet.h"
#include "angles_montecarlo.h"
#include "ccd_matr.h"
#include "telescope.h"
#include <thread>
#include <mutex>
#include <string>

using namespace std;

mutex mu_cout;
mutex mu_cout_thread_finished;

/*

// defined in a class Dust_cloud

#define ARRAY_X 80         ///< X dimensions for 3D array
#define ARRAY_Y 80         ///< Y dimensions for 3D array
#define ARRAY_Z 80         ///< Z dimensions for 3D array

#define star_quantity 10

struct star_coordinates
{
    float radius;       // radius of star
    float cloud_radius; // radius around star, there are no cosmic dust
    float x_pos;        // x - position of star
    float y_pos;        // y - positien of star
    float z_pos;        // z - position of star
};

*/


void photon_packages_traveling(int thread_number, Sub_Volume *arr_Vol[ARRAY_X][ARRAY_Y][ARRAY_Z], CCD_matr *matrix_obj,
                               star_coordinates *star_ar[star_quantity], short outer_Rad,  int photon_packet_q, unsigned int number_thre);

void cout_with_mutex(int thread_number, int pack_number);

void cout_thread_finished(int thread_number, short star_nr);


int main()
{
    short outer_Radius = 39; // sphere, when ARRAY_X = 80, ARRAY_Y = 80, ARRAY_Z = 80
    //short outer_Radius = 32; // assymetrical_1, when ARRAY_X = 80, ARRAY_Y = 80, ARRAY_Z = 80

    star_coordinates *star_arr[star_quantity];
    
    for(int i = 0; i < star_quantity; i++)
    {
        star_arr[i] = NULL; // because sometimes the value of unitialized pointer is undefined
    }
    
    star_arr[0] = new star_coordinates{1, 1, ARRAY_X/2, ARRAY_Y/2, ARRAY_Z/2};    // for sphere
    //star_arr[0] = new star_coordinates{1, 10, ARRAY_X/2, ARRAY_Y/2, ARRAY_Z/2};    // for sphere
    //star_arr[0] = new star_coordinates{1, 5, ARRAY_X/2 - 8, ARRAY_Y/2 - 8, ARRAY_Z/2}; // for assymetrical_1 (first star)
    //star_arr[1] = new star_coordinates{1, 5, ARRAY_X/2 + 8, ARRAY_Y/2 + 8, ARRAY_Z/2}; // for assymetrical_1 (second star)

    Dust_cloud *cloud_p;
    cloud_p = new Dust_cloud(1); // sphere
    //cloud_p = new Dust_cloud(2); // assymetrical_1
    cloud_p->set_min_max_phi_thet(0,0,2*PI,PI);
    cloud_p->set_outer_rad(outer_Radius);
    cloud_p->make_cloud_surface();

    cloud_p->cut_off_sphere(star_arr);
    cloud_p->fill_cloud();  

    
    Sub_Volume *vol_pointer_1, *vol_pointer_2, *vol_pointer_3, *vol_pointer_4, *vol_pointer_5, *vol_pointer_6,
               *vol_pointer_7, *vol_pointer_8, *vol_pointer_9, *vol_pointer_10, *vol_pointer_11, *vol_pointer_12;

    vol_pointer_1 = new Sub_Volume; // create only one object for all Sub_Volumes that is outer surface
    vol_pointer_1->set_cloud_surface(1);

    // for homogen cloud dust concentration = const
    //vol_pointer_1->set_scattering(1.7376); // scattering for homogen cloud
    //vol_pointer_1->set_extinction(1.871); // extinction for homogen cloud
    //vol_pointer_1->set_concentration(0.01843); // dust concentration for homogen cloud; TAU_0 = 1 optical way
    
    // for dust concentration / radius
    //vol_pointer_1->set_scattering(0.12854); // scattering for dust cloud, concentration / radius
    //vol_pointer_1->set_extinction(0.18363); // extinction for dust cloud, concentration / radius
    //vol_pointer_1->set_concentration(40); // dust concentration for cloud; TAU_0 = 10 optical way
    
    // for dust concentration / radius^2
    //vol_pointer_1->set_scattering(0.1343967); // scattering for dust cloud, concentration / radius^2
    //vol_pointer_1->set_extinction(0.1919953); // extinction for dust cloud, concentration / radius^2
    //vol_pointer_1->set_concentration(700); // dust concentration for cloud; TAU_0 = 10 optical way


    vol_pointer_2 = new Sub_Volume; // create only one object for all Sub_Volumes that volume between surfaces
    vol_pointer_2->set_cloud_surface(2);

    // for homogen cloud dust concentration = const
    //vol_pointer_2->set_scattering(1.7376); // scattering for homogen cloud
    //vol_pointer_2->set_extinction(1.871); // extinction for homogen cloud
    //vol_pointer_2->set_concentration(0.01843); // dust concentration for homogen cloud; TAU_0 = 1 optical way
    
    // for dust concentration / radius
    //vol_pointer_2->set_scattering(0.12854); // scattering for dust cloud, concentration / radius
    //vol_pointer_2->set_extinction(0.18363); // extinction for dust cloud, concentration / radius
    //vol_pointer_2->set_concentration(40); // dust concentration for cloud; TAU_0 = 10 optical way
    
    // for dust concentration / radius^2
    //vol_pointer_2->set_scattering(0.1343967); // scattering for dust cloud, concentration / radius^2
    //vol_pointer_2->set_extinction(0.1919953); // extinction for dust cloud, concentration / radius^2
    //vol_pointer_2->set_concentration(700); // dust concentration for cloud; TAU_0 = 10 optical way


    vol_pointer_3 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_3->set_cloud_surface(3);

    // for homogen cloud dust concentration = const
    //vol_pointer_3->set_scattering(1.7376); // scattering for homogen cloud
    //vol_pointer_3->set_extinction(1.871); // extinction for homogen cloud
    //vol_pointer_3->set_concentration(0.01843); // dust concentration for homogen cloud; TAU_0 = 1 optical way
    
    // for dust concentration / radius
    //vol_pointer_3->set_scattering(0.12854); // scattering for dust cloud, concentration / radius
    //vol_pointer_3->set_extinction(0.18363); // extinction for dust cloud, concentration / radius
    //vol_pointer_3->set_concentration(40); // dust concentration for cloud; TAU_0 = 10 optical way
    
    // for dust concentration / radius^2
    //vol_pointer_3->set_scattering(0.1343967); // scattering for dust cloud, concentration / radius^2
    //vol_pointer_3->set_extinction(0.1919953); // extinction for dust cloud, concentration / radius^2
    //vol_pointer_3->set_concentration(700); // dust concentration for cloud; TAU_0 = 10 optical way


    vol_pointer_4 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_4->set_cloud_surface(4);


    vol_pointer_5 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_5->set_cloud_surface(5);


    vol_pointer_6 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_6->set_cloud_surface(6);


    vol_pointer_7 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_7->set_cloud_surface(7);


    vol_pointer_8 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_8->set_cloud_surface(8);


    vol_pointer_9 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_9->set_cloud_surface(9);


    vol_pointer_10 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_10->set_cloud_surface(10);


    vol_pointer_11 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_11->set_cloud_surface(11);


    vol_pointer_12 = new Sub_Volume; // create only one object for all Sub_Volumes that is inner surface
    vol_pointer_12->set_cloud_surface(12);
    

    // array of pointers for sub_volume objects
    Sub_Volume *array_Vol[ARRAY_X][ARRAY_Y][ARRAY_Z];
    
    for(int i = 0; i < ARRAY_X; i++ )
    {
        for(int j = 0; j < ARRAY_Y; j++)
        {
            for(int k = 0; k < ARRAY_Z; k++)
            {
                array_Vol[i][j][k] = NULL;
            }
        }    
    }

    for(int i = 0; i < ARRAY_X; i++ )
    {
        for(int j = 0; j < ARRAY_Y; j++)
        {
            for(int k = 0; k < ARRAY_Z; k++)
            {
                switch(cloud_p->get_array_value(i,j,k))
                {
                    case 1:
                        array_Vol[i][j][k] = vol_pointer_1;
                        break;
                    case 2:
                        array_Vol[i][j][k] = vol_pointer_2;
                        break;
                    case 3:
                        array_Vol[i][j][k] = vol_pointer_3;
                        break;
                    case 4:
                        array_Vol[i][j][k] = vol_pointer_4;
                        break;
                    case 5:
                        array_Vol[i][j][k] = vol_pointer_5;
                        break;
                    case 6:
                        array_Vol[i][j][k] = vol_pointer_6;
                        break;
                    case 7:
                        array_Vol[i][j][k] = vol_pointer_7;
                        break;
                    case 8:
                        array_Vol[i][j][k] = vol_pointer_8;
                        break;
                    case 9:
                        array_Vol[i][j][k] = vol_pointer_9;
                        break;
                    case 10:
                        array_Vol[i][j][k] = vol_pointer_10;
                        break;
                    case 11:
                        array_Vol[i][j][k] = vol_pointer_11;
                        break;
                    case 12:
                        array_Vol[i][j][k] = vol_pointer_12;
                        break;
                }
            }
        }
    }

    delete cloud_p;

    CCD_matr *matr_pointer = new CCD_matr;
    
    int photon_pack_q = 1000000000; // how many photon packets will be scattered by programm

    //photon_packages_traveling( 0 ,array_Vol , matr_pointer, star_arr, outer_Radius,
    //                          photon_pack_q, 1 ); // only one thread, no multithreading


    //---- paralelizing, with several threads, choose the same thread number as number of processor cores ----//
    //
    thread thread_array[32];

    unsigned int number_threads = thread::hardware_concurrency(); // for multithreading, how many cores can be used

    if(number_threads > 32)
    {
        cout << "Core number is more than 32, only one thread will be created. \n";

        photon_packages_traveling( 0,array_Vol , matr_pointer, star_arr, outer_Radius,
                                      photon_pack_q, 1 ); // only one thread, no multithreading
    }
    else
    {
        cout << number_threads << " threads will be created \n";

        for(unsigned int numb_thread = 0; numb_thread < number_threads; numb_thread++ )
        {
            thread_array[numb_thread] = thread(photon_packages_traveling, 
                    numb_thread, array_Vol, matr_pointer, star_arr, outer_Radius, photon_pack_q, number_threads);
        }

        for(unsigned int numb_thread = 0; numb_thread < number_threads; numb_thread++)
        {
            thread_array[numb_thread].join();
        }
    }
    //

    matr_pointer->generate_gnu_txt_file();

    cout << "program works correctly \n";

    return 0;
}




void photon_packages_traveling(int thread_number, Sub_Volume *arr_Vol[ARRAY_X][ARRAY_Y][ARRAY_Z], CCD_matr *matrix_obj,
                               star_coordinates *star_ar[star_quantity], short outer_Rad, int photon_packet_q, unsigned int number_thre)
{
    short star_number = 0;
    
    Photon_source source = Photon_source();
    
    source.set_stokes_I_param(1);
    source.set_stokes_Q_param(0);
    source.set_stokes_U_param(0);
    source.set_stokes_V_param(0);

    //source.set_phot_weight(10000);
    source.set_phot_weight(1000); // for sphere

    Telescope telescope = Telescope();
    telescope.set_cloud_rad_AU(20); // set cloud radius as 20 AU
    telescope.set_dist_cloud_tel_PC(600); // set distance betwee cloud and telescope as 600 PC 
    telescope.set_telescope_rad_meters(4.1);    // set telescope radius as 4.1 meters
    telescope.calc_telescope_area();            // calculate telescope area
    telescope.calc_solid_angle();               // calculate telescope solid angle
    
    telescope.set_tel_polar_angle(PI/9*5);
    telescope.set_tel_azimuth_angle(PI/9*13);
    telescope.set_tel_dist_cloud_sys_tel_sys(2.4133e+8);
    telescope.calc_parameters_for_tel_param();

    //cout << "x = " << setprecision(15) << telescope.get_tel_X_coord() << "\n";
    //cout << "y = " << setprecision(15) << telescope.get_tel_Y_coord() << "\n";
    //cout << "z = " << setprecision(15) << telescope.get_tel_Z_coord() << "\n";

    Angles_montecarlo angles = Angles_montecarlo();

    emit_phot_param phot_structure;

    Photon_packet *packet;

    double max_optical_way = 0;
    double geom_way = 0;
    double optical_way = 0;
    int packet_quantity;
    
    //
    while(star_ar[star_number] != NULL)
    { 
        packet_quantity = photon_packet_q / number_thre; // how many photon packets will be scattered by one thread
        
        for(int i = 0; i < packet_quantity; i++)
        {
            if(i % 10000 == 0)
            {
                cout_with_mutex(thread_number, i);
            }

            source.initialize_photon_source(star_ar[star_number]);

            source.gener_emit_angles_point();  

            phot_structure = source.emit_phot();

            packet = new Photon_packet();
            
            /////////////// initializing photon packet with number of all stars ////////
            packet->initialize_struct_star_p(star_ar);
            packet->set_star_numb(star_number); // setting number of emitting star

            packet->initialize_packet(phot_structure); // setting PHI_1 , THET_1 and photon emission point
            packet->set_tel_parameters(telescope.get_tel_param()); // setting parameters of telescope
            packet->set_cloud_type(1); // 1 - setting cloud as sphere
                                   // 2 - setting cloud as assymetrical_1 
            packet->set_concentration_type(1); // 1 - concentration = const
                                                // 2 - concentration = concentration / radius
                                                // 3 - concentration = concentration / radius^2
            packet->set_CCD_matr_point(matrix_obj);
            packet->set_outer_inner_radius(outer_Rad, star_ar[star_number]->cloud_radius); 
                                                // essential to set before optical or geometrical way calculation
            max_optical_way = packet->max_optical_way_calc(arr_Vol);

            //cout << "max_optical_way = " << max_optical_way << "\n";

            optical_way = packet->gener_first_optical_way();

            //
            geom_way = packet->calc_geometrical_way(arr_Vol);

            if(packet->get_photon_location() != 3)
            {
                packet->calc_spher_coor_cur_tel(); // calculates spherical coordinates in direction to telescope
                packet->calc_Thet_peel();
                angles.set_angles_for_peel(packet->get_angles_for_peel());

                angles.set_iinc_qinc_uinc_vinc(1, packet->get_stokes_q_param(), packet->get_stokes_u_param(), packet->get_stokes_v_param());
                
                angles.calc_stokes_tel_par(); // calculates Stokes parameters for telescope: q_tel, u_tel, v_tel
                packet->set_stokes_q_tel(angles.get_q_tel());
                packet->set_stokes_u_tel(angles.get_u_tel());
                packet->set_stokes_v_tel(angles.get_v_tel());

                packet->set_F11_value( angles.get_F11_val() ); // sets F11 for Thet_peel value 
                packet->set_F12_value( angles.get_F12_val() ); // sets F12 for Thet_peel value

                packet->set_cos_2psi_tel(angles.calc_cos_2psi_tel());
                packet->set_sin_2psi_tel(angles.calc_sin_2psi_tel());

                packet->optical_way_peel_calc(arr_Vol); // calculate optical way for peeling off method                

                if(arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()][(int)packet->get_coord_Z_current()] != NULL)
                {
                    packet->first_scatter_phot_pack(arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()]
                                                        [(int)packet->get_coord_Z_current()]->get_scattering() /
                                                        arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()]
                                                        [(int)packet->get_coord_Z_current()]->get_extinction());
                }
                else
                {
                    packet->set_photon_location(3);
                }
                
            }

            //
            while(packet->get_photon_location() != 3 && packet->get_photon_location() != 4)
            {
                // angles calculation
                angles.calc_bigThet();
                angles.bisection(0.0001);
                angles.set_THET_1_PHI_1(packet->get_THET_1(), packet->get_PHI_1());
                packet->set_THET_1(angles.calc_THET_2());
                packet->set_PHI_1(angles.calc_PHI_2());
                angles.calc_stokes_scat_par();
                packet->set_stokes_q_param(angles.get_q_sca());
                packet->set_stokes_u_param(angles.get_u_sca());
                packet->set_stokes_v_param(angles.get_v_sca());
                optical_way = packet->gener_optical_way();

                geom_way = packet->calc_geometrical_way(arr_Vol);
                //packet->check_opt_way_simps(arr_Vol);
                //packet->check_opt_way_analytical(arr_Vol);

                if(packet->get_photon_location() != 3 )
                {  
                    packet->calc_spher_coor_cur_tel(); // calculates spherical coordinates in direction to telescope
                    packet->calc_Thet_peel();
                    angles.set_angles_for_peel(packet->get_angles_for_peel());
                    angles.set_iinc_qinc_uinc_vinc(1, packet->get_stokes_q_param(), packet->get_stokes_u_param(), packet->get_stokes_v_param());

                    angles.calc_stokes_tel_par(); // calculates Stokes parameters for telescope: q_tel, u_tel, v_tel
                    packet->set_stokes_q_tel(angles.get_q_tel());
                    packet->set_stokes_u_tel(angles.get_u_tel());
                    packet->set_stokes_v_tel(angles.get_v_tel());

                    packet->set_F11_value( angles.get_F11_val() ); // sets F11 for Thet_peel value
                    packet->set_F12_value( angles.get_F12_val() ); // sets F12 for Thet_peel value

                    packet->set_cos_2psi_tel(angles.calc_cos_2psi_tel());
                    packet->set_sin_2psi_tel(angles.calc_sin_2psi_tel());

                    packet->optical_way_peel_calc(arr_Vol); // calculate optical way for peeling off method                    

                    if(arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()][(int)packet->get_coord_Z_current()] != NULL)
                    {
                        packet->sca_phot_pack(arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()]
                                                        [(int)packet->get_coord_Z_current()]->get_scattering() /
                                                        arr_Vol[(int)packet->get_coord_X_current()][(int)packet->get_coord_Y_current()]
                                                        [(int)packet->get_coord_Z_current()]->get_extinction());
                    }
                    else
                    {
                        packet->set_photon_location(3);
                    }
                    
                }
            }
            //
        
            delete packet;
        }

        cout_thread_finished(thread_number, star_number);
        star_number++;
    }
    //
}


void cout_with_mutex(int thread_number, int pack_number)
{
    lock_guard<mutex> locker_2(mu_cout);

    cout << "Thread number " << thread_number << " = " << pack_number << "\n";
}


void cout_thread_finished(int thread_number, short star_nr)
{
    lock_guard<mutex> locker_3(mu_cout_thread_finished);
    
    cout << "Thread number " << thread_number << " for star number " << star_nr << " is finished.\n"; 
}
