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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    if(is_initialized)
        return;

//Number of particles
    num_particles = 100;
    default_random_engine gen;
//SD
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];

    normal_distribution<double> dist_x(x,std_x);
    normal_distribution<double> dist_y(y,std_y);
    normal_distribution<double> dist_theta(theta,std_theta);

    for(int i=0; i< num_particles ; i++) {
        Particle my_particle;
        my_particle.id=i;
        my_particle.x= dist_x(gen);
        my_particle.y=dist_y(gen);
        my_particle.theta =dist_theta(gen) ;
        my_particle.weight = 1.0;
        particles.push_back(my_particle);

    }
    is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

//According to formula applying the transformations
    for(int i=0; i< num_particles; i++)
    {
        if(yaw_rate < fabs(0.0001))
        {
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        }
        else
        {

            double yawDt = yaw_rate * delta_t;
            double outerTerm= velocity / yaw_rate;
            particles[i].x  = particles[i].x +  outerTerm * (sin(particles[i].theta + yawDt) - sin(particles[i].theta));
            particles[i].y  = particles[i].y +  outerTerm * (-cos(particles[i].theta + yawDt) + cos(particles[i].theta));
            particles[i].theta  = particles[i].theta + yawDt;

            particles[i].x = particles[i].x + dist_x(gen);
            particles[i].y = particles[i].y + dist_y(gen);
            particles[i].theta = particles[i].theta + dist_theta(gen);
        }
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    for ( int i = 0; i < observations.size(); i++) {
        LandmarkObs o = observations[i];
        double min_dist = 9999999;
        int min_id = -9999999;

        for (int j = 0; j < predicted.size(); j++) {
            LandmarkObs p = predicted[j];

            double e_dist = dist(o.x, o.y, p.x, p.y);

            if (e_dist < min_dist) {
                min_dist = e_dist;
                min_id = p.id;
            }
        }
        observations[i].id = min_id;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

//1. Transform observation from car coordinates to map coordinates
    //2. Euclidean function for landmarks and predicted particles KNN. In landmarks, considering the range of the sensor we need to calculate the distance of all particles int that range and chose the nearest one.
    //3. Distance between predicted and observed. Store all the ID's of the particles which are under minimum range defined by us.
    //4. Update weights as in Slide 19

//Step 1 Start
    vector<LandmarkObs> translated_landmarks;
    vector<LandmarkObs> predictedLandmarks;
    for(int i=0; i<num_particles; i++)
    {
        double t_x = cos(particles[i].theta)*observations[i].x - sin(particles[i].theta)*observations[i].y + particles[i].x;
        double t_y = sin(particles[i].theta)*observations[i].x + cos(particles[i].theta)*observations[i].y + particles[i].y;
        translated_landmarks.push_back(LandmarkObs{ observations[i].id, t_x, t_y });

        //}
// Step 1 End

//Step 2 Start

//for(int i=0;i<num_particles;i++)
//{
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double distance_x = map_landmarks.landmark_list[j].x_f - particles[i].x;
            double distance_y = map_landmarks.landmark_list[j].y_f - particles[i].y;
            if(fabs(distance_x)  <= sensor_range && fabs(distance_y) <=sensor_range)
            {

                predictedLandmarks.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
            }
        }

//Step 2 End

//Step 3 Start
        dataAssociation(predictedLandmarks, translated_landmarks);
// Step 3 End

//Step 4 Start
        double sig_x= std_landmark[0];
        double sig_y= std_landmark[1];
        particles[i].weight = 1.0;
        for(int k=0; k<translated_landmarks.size(); k++)
        {
            double x_obs= translated_landmarks[k].x;
            double y_obs= translated_landmarks[k].y;
            double mu_x;
            double mu_y;
            for(int j=0; j<predictedLandmarks.size(); j++)
            {
                if(predictedLandmarks[j].id == translated_landmarks[k].id)
                {
                    mu_x=predictedLandmarks[j].x;
                    mu_y=predictedLandmarks[j].y;

                }
            }

            double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
            double exponent= ( pow(mu_x-x_obs,2)/(2*pow(sig_x, 2)) + (pow(mu_y-y_obs,2)/(2*pow(sig_y, 2))) );
            double weightLocal = gauss_norm * exp(-exponent);
            particles[i].weight *= weightLocal;
        }


    }
//Step 4 End



}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> new_particles;
    vector<double> weights;
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    random_device rd;
    mt19937 gen(rd());

    discrete_distribution<int> uniintdist(0, num_particles-1);
    auto index = uniintdist(gen);

    // get max weight
    double max_weight = *max_element(weights.begin(), weights.end());

    // uniform random distribution [0.0, max_weight)
    uniform_real_distribution<double> unirealdist(0.0, max_weight); // check why discrete is not working here -- because it is not integral type it is double

    double beta = 0.0;

    // spin the resample wheel!
    for (int i = 0; i < num_particles; i++) {
        beta += unirealdist(gen) * 2.0;
        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        new_particles.push_back(particles[index]);
    }

    particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
        const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
