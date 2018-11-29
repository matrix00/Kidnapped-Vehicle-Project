/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

//have a random engine
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 50;
	is_initialized = true;


	//distribution for parameters
	normal_distribution<double> dist_x(x, std[0]);

	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);


	for (int i=0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.x= dist_x(gen);
		particle.y=dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight=1;
		
		particles.push_back(particle);	

		//cout<<"*************************int i=" << i << " x=" << particle.x << " y=" << particle.y << " theta=" << particle.theta << " weight " << particle.weight<<endl;
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	std::normal_distribution<double> xd(0,std_pos[0]);
	std::normal_distribution<double> yd(0,std_pos[1]);
	std::normal_distribution<double> theta_d(0,std_pos[2]);
	
	for (auto& p: particles) {

		//retrieve old values
		double x0 = p.x;
		double y0 = p.y;
		double theta_0 = p.theta;

		if (fabs(yaw_rate) < 0.0001) {

			//predict new values
			 p.x = x0+velocity * delta_t * cos(theta_0);
			 p.y = y0+velocity * delta_t * sin(theta_0);
		}
		else {
			//predict new values
			 p.x = x0+velocity * (sin(theta_0+yaw_rate*delta_t)-sin(theta_0))/yaw_rate;
			 p.y = y0+velocity * (cos(theta_0) - cos(yaw_rate*delta_t+theta_0))/yaw_rate;
			 p.theta = theta_0+yaw_rate*delta_t;
		}
	
		//generate random noise for each particle
		double noise_x = xd(gen);
		double noise_y = yd(gen);
		double noise_theta = theta_d(gen);

		//add noise to prediction
		p.x  += noise_x;	
		p.y  += noise_y;	
		p.theta += noise_theta;	

		//cout<<"predict " << " x=" << p.x << " y=" << p.y << " theta=" << p.theta << " weight " << p.weight<<endl;
	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//cout<<"dassoc predicted=" << predicted.size() << " obs size"<< observations.size()<<endl;

	
	for (auto& ob:observations){
		double min = numeric_limits<double>::max();
		int min_id = -1;
		//find nearest neighbor
		for (auto& pred: predicted){
			double distance = dist(pred.x, pred.y, ob.x, ob.y);
			if (min > distance) {
				min = distance;
				min_id = pred.id;
			}
		}
		ob.id = min_id;
	}
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for (auto& p: particles) {

		//predict the landmark with sensor range

		vector<LandmarkObs> predictions;
		for (auto& lmark:map_landmarks.landmark_list) {
			double pldist = dist(p.x, p.y, lmark.x_f, lmark.y_f);
			if (pldist < sensor_range)
				predictions.push_back(LandmarkObs{lmark.id_i, lmark.x_f, lmark.y_f});
		}


        	//transformation requires both rotation AND translation (but no scaling).

		vector<LandmarkObs> xformedObs;
		for (auto& obs:observations) {
			//http://planning.cs.uiuc.edu/node99.html
			LandmarkObs xformedOb;
			xformedOb.x = obs.x* cos(p.theta) - obs.y * sin(p.theta) + p.x;
			xformedOb.y = obs.x* sin(p.theta) + obs.y * cos(p.theta) + p.y;
			xformedOb.id = obs.id;	
			xformedObs.push_back(xformedOb);
		}
		
		// use data association to find closed one

		dataAssociation(predictions, xformedObs);

		p.weight=1.0;

		//update weight
		for (auto& xObs: xformedObs) {
			//get the predicted one that matches observation
			LandmarkObs pmatch;			
			for (auto& pred: predictions) {
				if (pred.id == xObs.id) {
					pmatch = pred;
					//cout <<"########################### found ###############"<< pmatch.id << endl; 
				}

			}

			//https://en.wikipedia.org/wiki/Multivariate_normal_distribution
			//using multi variate gaussian distribution

			// f(x,y) = 1/2 PI SigmaX SigmaY exp(- 1/2 (x - mu X)**2/ SigmaX**2 + (y - mu Y)**2/ SigmaY** 2)

			double obsWeight = (1.0/(2*M_PI * std_landmark[0] * std_landmark[1])) * exp(-0.5*((pow(xObs.x-pmatch.x, 2)/pow(std_landmark[0], 2))+ (pow(xObs.y-pmatch.y, 2)/pow(std_landmark[1],2))));
			
			p.weight *= obsWeight;
			//cout<<"Obs weight " << obsWeight <<  " pweight " << p.weight <<endl;
		
		}

		weights.push_back(p.weight);
	
	}
	
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//use discrete distribution for resampling
	//reference is helpful. Much easier than running the wheel for resampling
	vector<Particle> resampledPs;

	std::random_device rd;
	std::mt19937 gen(rd());
    	std::discrete_distribution<> d(weights.begin(), weights.end());

	for (int i=0; i < num_particles; i++) {
		int index = d(gen);
		resampledPs.push_back(particles[index]);
	}

	particles = resampledPs;

	//cout << "resampled " << particles.size() << endl;
	weights.clear();
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
