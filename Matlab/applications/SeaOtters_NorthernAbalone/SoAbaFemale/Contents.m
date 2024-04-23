% SOABAFEMALE
%
% Files
%   Abalone_study                 - Abalone population study
%   Abalone_study_adults          - Abalone population study
%   abaloneDensity2state          - function abaloneState=abaloneDensity2state(d_aba)
%   build_joint_reward            - function R=build_joint_reward(k_aba,level_so,level_aba)
%   build_state_matrix            - function state_matrix= build_state_matrix(nbs_aba,nbs_seaotter)
%   compute_expert_policy         - function expi=compute_expert_policy(expert_num)
%   compute_poaching_impact       - function N=compute_poaching_impact(N,poach)
%   compute_transition            - function Tr=ComputeTransition(iFR,iFRnum,ksim)
%   derive_hyp_FR                 - function [cmin,cmax,dmin,dmax]=derive_hyp_FR(k_aba_max)
%   derive_sig_FR                 - function [cmin,cmax,dmin,dmax]=derive_sig_FR(k_aba_max)
%   draw_policy                   - function draw_policy(Pi_Best,fname,liste)
%   evaluate_strategy             - function evaluate_strategy(policyPi,str,iFR,iFRnum)
%   evaluate_strategy_mean        - function rewardSum=evaluate_strategy_mean(PolicyPi,str,iFR,iFRnum,nsim)
%   explore_Q                     - function liste_a=explore_Q(PiB,sensitivity)
%   figure_s2                     - 
%   get_results                   - function All=get_results(name)
%   initialising_northern_abalone - function N1=initialising_northern_abalone(disp,target_aba_density)
%   load_param                    - function load_param( isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%   main_experts                  - function main_experts(nsim,               % simulation number
%   main_SDP                      - function main_SDP(nsim,                % number of simulation
%   main_SDP_SA                   - function main_SDP_SA(nsim,                % number of simulation
%   northern_abalone_growth_t     - function Ni1=northern_abalone_growth_t(Ni1)
%   perform_action                - function [cost,nbso,oil,AbaPop]=
%   plot_joint_reward_fig1        - function R=build_joint_reward(k_aba,level_so,level_aba)
%   predation_FR                  - Function new_Tabundance_prey=
%   results_analysis              - function results_analysis(filedate,EOA)
%   results_analysis_SA           - function results_analysis_SA(filedate,EOA)
%   sea_otter_growth              - Function [Y,OS]=sea_otter_growth(N)
%   sea_otters_study              - Sea otter population model including risk of oil spill
%   seeIndex                      - function index= seeIndex(state)
%   seeState                      - function state= seeState(index,state_matrix)
%   sensitivity                   - sensitivity analysis script
%   show_all_FR                   - 
%   show_FR                       - function function show_FR(iFR,iFRnum)
%   simulation_t                  - function
%   SOabundance2state             - Function seaotterState = SOabundance2state(popSeaotter)
%   table_supp_info               - for a given k and r provides table
%   table_supp_info2              - provides table s3 - use results_analysis_SA
%   value_iteration               - function Pib=value_iteration()
