/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Survey.hpp
 * Author: matthewsupernaw
 *
 * Created on January 30, 2017, 3:22 PM
 */

#ifndef SURVEY_HPP
#define SURVEY_HPP

#include "Selectivity.hpp"

namespace mas{
    
    template<typename REAL_T>
    struct Survey : public mas::ModelObject<REAL_T>{
        int id;
        std::string name;
        int population;
        int selectivity_model_id;
        std::vector<int> area_ids;
        std::shared_ptr<mas::SelectivityBase<REAL_T> > selectivity;
         std::unordered_map<int, std::unordered_map<int, int> > area_season_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > area_season_selectivity;

        std::unordered_map<int, std::unordered_map<int, int> > season_area_selectivity_ids;
        std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > > season_area_selectivity;
        
         typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_selectivity_ids_iterator;
        typedef typename std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > >::iterator area_sectivity_iterator;
        typedef typename std::unordered_map<int, std::unordered_map<int, int> >::iterator season_area_id_iterator;
        typedef typename std::unordered_map<int, int>::iterator area_id_iteraor;
        typedef typename std::unordered_map<int, int>::iterator season_id_iteraor;
        typedef typename std::unordered_map<int, std::unordered_map<int, std::shared_ptr<mas::SelectivityBase<REAL_T> > > >::iterator season_area_selectivity_iterator;
        
    };
    
    template<typename REAL_T>
    std::ostream& operator << (std::ostream& out, const Survey<REAL_T>& survey){
        out <<"Survey\n";
        out<<"Name: " <<survey.name<<"\n";
        out<<"Areas: [ ";
        
        out<<"# Selectivity: " <<survey.area_season_selectivity_ids.size()<<"\n";
        return out;
                
    }
    
}


#endif /* SURVEY_HPP */

