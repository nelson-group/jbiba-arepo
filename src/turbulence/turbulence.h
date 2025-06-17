#ifdef TURBULENCE
double turb_random_normal(double mean, double variance);
void turb_OU_check_update();
void OU_update();
void get_driving_vec(double pos[3], double acc3d[3]);
void turb_apply_driving();
#endif