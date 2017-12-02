# **Unscented Kalman Filter Project**

## **1.compile**
The code compiled without errors. The result is in "[build](./build)" folder.

## **2.Accuracy**
The RMSE = [0.0848, 0.0889, 0.3720, 0.2992]  
[Here](./video/ukf.ogv) (./video/ukf.ogv)is the vido of my monitor.

To get the accuracy I tried a lot of combinations of the `std_a_` and `std_yawdd_`, considering the NIS and the RMSE.

Instead of ploting the NIS, I calculated the ratio of the NIS over 5% value.

Finally I got the value `std_a_=3` and `std_yawdd_=0.2`.

The ratio of my NIS is shown below (it looks like reasonable):

|senser|ratio|
|------|-----|
|laser|10/249(0.0401606)|
|radar|10/249(0.0401606)|


## **3.The follows of my algorithm**
### *3.1 I use the first measurements to initialize the `x_` and `P_`.*

According to the sensor type, I used diffrent values to initialize the `x_` and `P_`. 

For laser:  
```c++
x_<<meas_package.raw_measurements_(0),
  meas_package.raw_measurements_(1),
  0,
  0,
  0;
P_<<std_laspx_*std_laspx_,0,0,0,0,
    0,std_laspy_*std_laspy_,0,0,0,
    0,0,1,0,0,
    0,0,0,1,0,
    0,0,0,0,1;
```
For radar:
```c++
x_<<meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1)),
meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1)),
0,
meas_package.raw_measurements_(1),
0;
P_<<std_radr_*std_radr_,0,0,0,0,
    0,std_radr_*std_radr_,0,0,0,
    0,0,std_radrd_*std_radrd_,0,0,
    0,0,0,std_radphi_*std_radphi_,0,
    0,0,0,0,1;
```

### *3.2 After the first measurement, I make predicion and update step depends on the sensor type.*


```c++
double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    Prediction(dt);
    UpdateLidar(meas_package);
    time_us_=meas_package.timestamp_;
  }
  if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    Prediction(dt);
    UpdateRadar(meas_package);
    time_us_=meas_package.timestamp_;
  }
```

## **4.Catch the car**

The catch project is in the folder named "Bonus"

[Here](./video/catch.ogv)(./video/catch.ogv) is the vido of my result.

