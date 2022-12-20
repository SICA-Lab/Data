% this code work for evalatue TA source amplitude
% for different saturation levels
%% input parameters, 1 = water, 2 = rock
clear;close all; clc;
eps_1 = 80;
sig_1 =0.05;
rho_1 = 1000;
alpha_1 = 210e-6;
C_p1 = 4180;
C_l1 = 1500;
C_s1 = 0;
K_1 = C_l1^2*rho_1;

eps_2 = 4.27;
sig_2 = 0;
rho_2 = 2650;
alpha_2 = 1e-6;
C_p2 = 743;
C_l2 = 6008;
C_s2 = 4075;
mu_2 = rho_2*C_s2^2;
K_2 = C_l2^2*rho_2 - 4/3*mu_2;

Phi = [0.06:0.06:0.3];
f = 1.3e9; % EM frequency
my_eps = zeros(size(Phi))';
my_sig = zeros(size(Phi))';
my_Cp = zeros(size(Phi))';
my_alpha = zeros(size(Phi))';
Q = zeros(size(Phi))';
Q1 = zeros(size(Phi))';

for ii = 1:length(Phi)
    phi = Phi(ii);% porosity
    eps_c1 = eps_1 - j*sig_1/(2*pi*f);
    eps_c2 = eps_2 - j*sig_2/(2*pi*f);
    syms eps_c
    eqn =  ((eps_c - eps_c2)/(eps_c1 - eps_c2))*(eps_c1/eps_c)^(1/3) == phi;
    S = vpasolve(eqn,eps_c);
    
    alpha = ((1-phi)*K_2*alpha_2 + phi*K_1*alpha_1)/((1-phi)*K_2 + phi*K_1);
    C_p = (rho_2*C_p2*(1-phi) + rho_1*C_p1*phi)/rho_2;
    
    my_eps(ii) = real(S);
    my_sig(ii) = -imag(S)*2*pi*f;
    my_Cp(ii) = C_p;
    my_alpha(ii) = alpha;
end
SRC = [];
%% computation

for ii = 1:length(Phi)
    clear model;
    import com.comsol.model.*
    import com.comsol.model.util.*
    model = ModelUtil.create('Model');
    model.modelPath('C:\Users\Chang\Desktop');
    model.label('test1.mph');
    model.hist.disable;
    
    model.component.create('comp1', true);
    
    model.component('comp1').geom.create('geom1', 3);
    
    model.result.table.create('tbl1', 'Table');
    model.result.table.create('tbl2', 'Table');
    
    model.component('comp1').mesh.create('mesh1');
    
    
    model.component('comp1').geom('geom1').lengthUnit('mm');
    model.component('comp1').geom('geom1').create('blk1', 'Block');
    model.component('comp1').geom('geom1').feature('blk1').set('size', [165 82.5 100]);
    model.component('comp1').geom('geom1').create('blk4', 'Block');
    model.component('comp1').geom('geom1').feature('blk4').set('pos', [0 0 100]);
    model.component('comp1').geom('geom1').feature('blk4').set('size', [165 82.5 6]);
    model.component('comp1').geom('geom1').create('blk2', 'Block');
    model.component('comp1').geom('geom1').feature('blk2').set('pos', [0 0 100+6]);
    model.component('comp1').geom('geom1').feature('blk2').set('size', [165 82.5 80]);
    model.component('comp1').geom('geom1').create('blk3', 'Block');
    model.component('comp1').geom('geom1').feature('blk3').set('pos', [67.5 34.75 125+6]);
    model.component('comp1').geom('geom1').feature('blk3').set('size', [30 13 30]);
    model.component('comp1').geom('geom1').run;
    
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material.create('mat2', 'Common');
    model.component('comp1').material.create('mat3', 'Common');
    model.component('comp1').material.create('mat4', 'Common');
    model.component('comp1').material('mat2').selection.set([2 3]);
    model.component('comp1').material('mat3').selection.set([4]);
    model.component('comp1').material('mat4').selection.set([2]);
    
    model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
    model.component('comp1').physics('emw').create('port1', 'Port', 2);
    model.component('comp1').physics('emw').feature('port1').selection.set([3]);
    model.component('comp1').physics('emw').create('sctr1', 'Scattering', 2);
    model.component('comp1').physics('emw').feature('sctr1').selection.set([4 5 7 8 10 12 13 21 22]);
    model.component('comp1').physics('emw').feature('port1').set('Pin', '1000 [W]');
    model.component('comp1').physics('emw').feature('port1').set('PortType', 'Rectangular');
    
    model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
    model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1');
    
    model.result.table('tbl2').comments('Volume Average 1 (emw.Qh)');
    
    model.component('comp1').material('mat1').label('Air');
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
    model.component('comp1').material('mat2').label('Oil');
    model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'2.6' '0' '0' '0' '2.6' '0' '0' '0' '2.6'});
    model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0.01' '0' '0' '0' '0.01' '0' '0' '0' '0.01'});
    model.component('comp1').material('mat3').label('Sand');
    model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity',my_sig(ii));
    model.component('comp1').material('mat3').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', my_eps(ii));
    model.component('comp1').material('mat4').label('Acrylic ');
    model.component('comp1').material('mat4').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'3.4' '0' '0' '0' '3.4' '0' '0' '0' '3.4'});
    model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'6.7e-15' '0' '0' '0' '6.7e-15' '0' '0' '0' '6.7e-15'});
    
    
    
    model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
    model.component('comp1').mesh('mesh1').feature('size').set('hmax', 5);
    model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.5);
    model.component('comp1').mesh('mesh1').run;
    
    model.study.create('std1');
    model.study('std1').create('freq', 'Frequency');
    
    
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('p1', 'Parametric');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    
    
    model.study('std1').feature('freq').set('plist', '1.3[GHz]');
    
    model.sol('sol1').attach('std1');
    model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
    model.sol('sol1').feature('v1').set('cname', {'freq'});
    model.sol('sol1').feature('v1').set('clist', {'1.3[GHz]'});
    model.sol('sol1').feature('s1').set('stol', 0.01);
    model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
    model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
    model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'1.3[GHz]'});
    model.sol('sol1').feature('s1').feature('p1').set('punit', {'GHz'});
    model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
    model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
    model.sol('sol1').feature('s1').feature('i1').label('Suggested Iterative Solver (emw)');
    model.sol('sol1').feature('s1').feature('i1').set('itrestart', 300);
    model.sol('sol1').feature('s1').feature('i1').set('prefuntype', 'right');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('iter', 1);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('sorvecdof', {'comp1_E'});
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('sorvecdof', {'comp1_E'});
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
    model.sol('sol1').runAll;
    
    model.result.numerical.create('av1', 'AvVolume');
    model.result.numerical('av1').selection.set([4]);
    model.result.numerical('av1').set('probetag', 'none');
    model.result.numerical('av1').set('table', 'tbl2');
    model.result.numerical('av1').set('expr', {'emw.Qh'});
    model.result.numerical('av1').set('unit', {'W/m^3'});
    model.result.numerical('av1').set('descr', {'Total power dissipation density'});
    model.result.numerical('av1').setResult;
    Q(ii) = model.result.numerical('av1').getReal;
    Q1(ii) = my_alpha(ii)/my_Cp(ii)*Q(ii);
 
    model.result.dataset.create('cpt1', 'CutPoint3D');
    model.result.create('pg1', 'PlotGroup1D');
    model.result('pg1').create('ptgr1', 'PointGraph');
    model.result.export.create('plot1', 'Plot');
    model.result.dataset('cpt1').label('Meas');
    model.result.dataset('cpt1').set('method', 'grid');
    model.result.dataset('cpt1').set('gridx', 'range(67.5,1,97.5)');
    model.result.dataset('cpt1').set('gridy', 'range(34.75,1,47.75)');
    model.result.dataset('cpt1').set('gridz', 'range(131,1,161)');
    model.result('pg1').label('E');
    model.result('pg1').set('data', 'cpt1');
    model.result('pg1').set('xlabel', 'freq (GHz)');
    model.result('pg1').set('ylabel', 'Electric field norm (V/m)');
    model.result('pg1').set('xlabelactive', false);
    model.result('pg1').set('ylabelactive', false);
    model.result.export('plot1').set('filename', strcat('C:\Users\Chang\Desktop\E',num2str(ii),'.txt'));
    model.result.export('plot1').set('header', false);
    
     t = mphinterp(model,'emw.Qh','dataset','cpt1');
     t = my_alpha(ii)/my_Cp(ii)*t;
     SRC = [SRC;t];

    mphsave(model,strcat('my8_',num2str(ii),'.mph'));
end
%% output
Q3 = [Q1];
save('Q3.mat','Q3');
save('SRC.mat','SRC');
