--------------------------------------------------------------------------
Copyright (c) 2011, CSIRO
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
README.txt file for the sea otter and abalone conservation management problem

RELEASE
    27/10/2020 Added to github
    15/05/2012 AbaSO Second release v. 0.2 
        - Changed Pemax value and fixed northern abalone population model.
    22/06/2011 AbaSO Initial release v. 0.1
    
INSTALLATION
    a- Extract the files from the archive file using winrar or winzip. The 
    files will be extracted into a directory AbaSO
    b- Open matlab, go to the directory AbaSO with all the matlab files.

GETTING STARTED
    Please refer to paper Chades et al (2012) for modelling details. 
	https://conbio.onlinelibrary.wiley.com/doi/abs/10.1111/j.1523-1739.2012.01951.x
    
    * If you wish to reproduce the figures from Chades et al (2011), try:
    >> results_analysis('24-Jun-2011','A')
    This function should give 3 graphs representing performances of the 
    optimal and experts strategies.
    
    * If you wish to derive performance results for different carrying 
    capacity parameters or growth rate parameters feel free to use:

    >> main_experts()

    or main_experts(nsim,               % number of simulations to evaluate strategies
                    isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
                    isDisplayOn,         % 0 or 1    figures off/on   (default 0)
                    maxKAba,             % carrying capacity in best habitat (default 3.34)
                    maxrAba,             % max growth rate in best habitat (default 1.6)
                    KSO,                 % carrying capacity sea otter (default 4073)
                    rSO)                 % growth rate sea otter (default 0.191)
    main_experts is the main function to simulate the performances of
    the predifined strategies. main_expert simulates and evaluates the 
    predefined strategies (4) under different functional responses 
    (linear, sigmoid, hyperbolic).
     
    EXAMPLE: >> main_experts(20,1,1,3.3,1.7,4073,0.191)

    * If you wish to compute a near optimal policy try:
    >> main_SDP()
    main_SDP is the main function to run the optimisation procedure.
    with main_SDP( nsim,                % number of simulations to learn the transition probabilities for each state action pair (default 500)
                   isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
                   isDisplayOn,         % 0 or 1    figures off/on   (default 0)
                   maxKAba,             % carrying capacity in best habitat (default 3.34)
                   maxrAba,             % max growth rate in best habitat (default 1.6)
                   KSO,                 % carrying capacity sea otter (default 4073)
                   rSO)                 % growth rate sea otter (default 0.191)
     
    For more information try:
    >> help main_SDP

    * Finally you can change most parameters used in this program by 
    editing the file load_param.m, including the functional responses.
with load_param( isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
                 isDisplayOn,         % 0 or 1    figures off/on   (default 0)
                 maxKAba,             % carrying capacity in best habitat (default 3.34)
                 maxrAba,             % max growth rate in best habitat (default 1.6)
                 KSO,                 % carrying capacity sea otter (default 4073)
                 rSO)                 % growth rate sea otter (default 0.191)
  
CITING/REFERENCE
    If you wish to refer to this work please cite:
    Chades I, Martin TG, Curtis J. (2012) Setting realistic recovery targets for interacting endangered species
    Conservation Biology. Or go to http://www.csiro.au/people/Iadine.Chades.html

KNOWN BUGS
   

TO DO
    Provide a GUI to set up parameters;
        
QUESTIONS/COMMENTS/CONCERNS

  If you have any questions, suggestions, problems or complaints
  regarding AbaSO, please email iadine chades <iadine.chades@csiro.au>
-------------------------------------------------------------------------------
end of README.txt