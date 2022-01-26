clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

%%

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
selectedpoint_str = sprintf('data/bone/amode_measure3.mat');
load(selectedpoint_str);
U = [vertcat(amode_prereg.Position); vertcat(amode_mid.Position)]';

% obtain all combination of z rotation and translation
range = 8;
step  = 0.25;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to translation vector
ts = [ zeros(2, length(t_z)); t_z];

% these are noise that is so difficult 
random_noise_special = [];
random_noise_special(:,:,1) = [ ... 
                                -0.000898170734688013	0.000261809513504595	0.000370534080999302	0.000321210334399039	0.000868345425639353	0.000480635772951166	7.63088439796986e-05	0.000870649765431192	-0.000796475731238937	-0.000467804997314066	0.000786014743272384	-0.000368788922298350	0.000929095782167068	-0.000331914890856042	0.000436758056871135; ...
                                -0.000887808513645899	0.000783955630026238	0.000391489586394435	3.98731907032893e-05	-0.000507046103701004	-0.000491662934057967	0.000820272268125519	0.000855438918265017	-0.000266248120555979	-0.000546922701145471	-9.30049855086274e-05	0.000988033750132379	0.000332530478104338	4.55455327057414e-05	0.000556035112715205; ...
                                -0.000329605200162935	0.000346793252868116	0.000599662866307252	-0.000336527178398671	2.24127422847769e-05	0.000691169517169519	-0.000298587179707060	0.000175400720791904	-0.000448308930064015	-0.000997358439853377	0.000156842011008397	0.000967785767106008	0.000452316352749056	-0.000453291268783303	-0.000837871560941901; ...
                                ] ;
random_noise_special(:,:,2) = [ ...
                                -0.000734033247567636	0.000863405551349943	-0.000707875719304627	-0.000912553202260724	-0.000126510361874170	-0.000549594181248037	-0.000494963191055261	0.000881012180201516	0.000302099285962251	-0.000241752538419689	0.000298545477653930	0.000626878259333601	-0.000685224904499018	0.000315407810100263	-0.000742780852935327; ...
                                -0.000859385518010654	-0.000782198172062283	0.000387925857128569	0.000595847653346224	-0.000310502291308271	0.000162025837460293	0.000308683319174011	-0.000814938844465143	7.82155272671518e-06	0.000537814825677038	0.000844367644759472	0.000575271190776897	-0.000170003456826103	0.000359954993533375	0.000527812762659028; ...
                                0.000412162883472582	0.000303393683538336	0.000251227245307360	0.000972421398167656	-0.000902387763151955	0.000202750676528279	-0.000879808808728803	-0.000332044137646533	6.18502812568969e-05	-0.000434977468440233	-0.000338656228760807	-0.000751355041111906	0.000374665596946727	-0.000981226645981506	0.000683830995946568; ...
                                ] ;
random_noise_special(:,:,3) = [ ...
                                -0.000402289131928562	0.000436919050210164	-6.35882060939781e-05	0.000824090888558019	0.000105262002129380	0.000203755431879534	0.000396714120476944	-0.000156428365873938	0.000661457641406877	0.000875215033780837	0.000621717062498548	0.000450978110814705	-0.000506333085188846	-9.38551633204507e-05	-0.000983643395038746; ...
                                0.000427361883789081	0.000990647466055207	-0.000812393884849246	-0.000827171133080045	0.000370298075305909	-0.000602992189521661	-8.06392391973229e-05	0.000192108732992176	-0.000754629983758925	0.000310416418030925	-0.000903146563627624	-0.000721911977928790	0.000679944546441567	-0.000210482527346529	-0.000203122620717038; ...
                                -0.000279027599164678	-0.000832599790018667	0.000450340912275219	3.99017774101302e-05	-0.000400690433750867	0.000313514830444218	-0.000684681661104178	-0.000355337007376920	-0.000497638596073755	0.000505821738043894	-0.000170581724605588	0.000265114518304006	0.000673341646582557	0.000923190603308830	-2.44554521771054e-05; ...
                                ] ;
random_noise_special(:,:,4) = [ ...
                                -0.000941846360469744	-0.000600719438227788	-0.000252806735538165	0.000403891378420410	-0.000500452529632527	-0.000842730650273481	-0.000953557604806731	0.000311896339067111	0.000791895559244191	-0.000171279823921665	-0.000279113013475847	0.000663699761329612	-0.000469216770349287	0.000136517176685086	0.000975903935312048; ...
                                0.000350780761869956	6.03230772211053e-05	0.000212907268038431	-0.000214252171159470	0.000454693466222801	-0.000395538143451128	-0.000647771496684588	0.000713075202643242	-0.000582419208255237	-0.000540869572272876	3.04813629886484e-06	-0.000211058035917689	0.000583890025104214	0.000876043114464499	0.000524697005050533; ...
                                -0.000776819253233500	0.000254713694460240	0.000692777369086583	0.000143737326051646	-0.000313341709281208	-0.000399195036152005	0.000721837593222053	-0.000458119068375827	0.000928336360492117	0.000636740397929515	0.000495348236872426	-0.000806727306236569	0.000873721485200616	0.000194496789100008	-0.000628768262786288; ...
                                ] ;
random_noise_special(:,:,5) = [ ...
                                0.000965375605948196	-7.57886910024864e-06	0.000922438024821965	0.000342533298063960	-0.000883584351006005	-0.000659060758358385	-0.000296025356686306	-6.31351419300616e-05	-0.000157128228260691	5.02920175326996e-05	0.000496660525606954	0.000372935467808939	-0.000315398274901000	0.000267092810234265	-0.000157082602281643; ...
                                0.000364442323193106	-0.000424675578742543	0.000512522951419713	-0.000627890779749577	-0.000911867555179792	0.000354201119539353	-0.000799682967835711	0.000374682413787116	0.000318685052241927	-0.000617437100355694	-0.000160346166667934	0.000911452025085338	-0.000225932855988662	0.000130502118090414	-7.55994559413815e-05; ...
                                0.000741801436129655	-0.000633961365225665	0.000900706981027967	-0.000420316083713807	-0.000440882337349986	0.000804835856652461	0.000586004471468939	0.000366157057494832	0.000266854258418611	2.38120194816587e-05	0.000474958838955568	-0.000384187077339670	-0.000947648189283128	-0.000268977174337537	-0.000371122447885548; ...
                                ] ;
random_noise_special(:,:,6) = [ ...
                                0.000540174904091290	-0.000488554546946192	-0.000318097511897920	-0.000559997834332685	0.000276433538768396	0.000536904227916605	-0.000895117530945080	0.000260459534513018	0.000233352463799721	0.000753659114487638	-0.000663859795377970	0.000306979709661502	-0.000401971592168668	-0.000938287360242694	-0.000867513416398698; ...
                                0.000882247934857022	-0.000243378384566586	0.000799104203005330	0.000982152799772396	8.26168229635339e-06	0.000568725775643616	-0.000353162032049570	0.000965595295907012	-0.000382565564686187	0.000181402877006034	-0.000320246493822583	-0.000511428851633316	-0.000162369810750932	-0.000626839231907595	-0.000333012615074606; ...
                                -0.000737070424431286	0.000987471391121623	-0.000525036775795841	0.000902210104162288	-0.000283763307583509	-0.000942163814915933	0.000599373249609289	-0.000683500526403846	-0.000829842987177013	-0.000393753725167676	-0.000865107772416430	0.000515836466493554	-0.000888335612490001	-0.000458188701559318	-0.000900188854588937; ...
                                ] ;
random_noise_special(:,:,7) = [ ...
                                -0.000720357956384741	-0.000951228696738581	0.000538657007158497	-3.73547526116666e-05	0.000157556779475063	-0.000715667632907684	-0.000744464329444152	-0.000872054963491580	0.000709783895756316	0.000761208585023159	-0.000312072623436067	0.000626476915375111	0.000259003343121730	0.000206794791004276	0.000475861612861615; ...
                                0.000621707569197752	6.14471631695316e-05	-0.000190495375357003	-0.000117648414464086	0.000567219400217990	0.000385837069673255	0.000281808794863282	4.70743462368081e-05	-0.000981969780506743	-0.000563900598912944	0.000600154794346928	-0.000827640789279883	-0.000872273090497455	0.000596888509533036	0.000514079036951721; ...
                                -0.000494280654608258	0.000869709348117818	0.000435621850572869	0.000884721657059260	0.000429614100948209	-0.000212702094191567	-4.98303069093207e-06	-0.000370920302857643	-0.000594829039023257	-0.000364620631647274	0.000142168926110751	-0.000615725026120594	0.000485924631331904	-0.000678810493333817	0.000762622184704587; ...
                                ];
random_noise_special(:,:,8) = [ ...
                                0.000477535427617170	-0.000892611765619094	0.000968089475979994	-0.000135992818150055	-0.000871803989438377	-0.000705663004808595	-0.000824211607027359	-0.000887400277446113	-0.000304553614608315	-0.000106242124142584	0.000875644237431087	0.000985591164900366	0.000732993171210047	-0.000981952542558144	0.000803244377612372; ...
                                -0.000680254482290720	0.000428779054167571	0.000186196392532406	0.000861024684685878	0.000318222625613707	0.000550456179813374	0.000800921391902665	0.000280112506752310	0.000333439012092989	-0.000965937953655870	-4.24518260419149e-05	-0.000191730748099891	0.000378046138154629	-0.000661845744923359	0.000418279794833208; ...
                                0.000462941430967641	-0.000756528715401839	-0.000460913589270358	1.57706170805919e-05	-0.000101655964310372	6.25235136876650e-05	0.000736302368051879	0.000778529259588922	-9.10785142166128e-05	0.000722835312392171	0.000139662987586339	-0.000560317362494831	-0.000657061233407718	-0.000587224458096992	-0.000736139900544583; ...
                                ] ;

noise        = 1;
num_trials   = 1;
costfunctions_min = zeros(num_trials, 2); 

for trial=1:num_trials
    
    fprintf('trials: %d\n', trial);

    % add isotropic zero-mean gaussian noise to U, simulating noise measuremen
    % random_noise = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
    random_noise = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(3,size(U, 2));
    model_ptCloud = (U + random_noise)';
%     model_ptCloud = (U + random_noise_special(:,:,1))';

%     % plot Ŭ, the noiseless, complete, moving dataset
%     figure1 = figure(1);
%     figure1.WindowState  = 'maximized';
%     axes1 = axes('Parent', figure1);
%     plot3( axes1, ...
%            U_breve(1,:), ...
%            U_breve(2,:), ...
%            U_breve(3,:), ...
%            '.', 'Color', [0.7 0.7 0.7],...
%            'MarkerSize', 0.1, ...
%            'Tag', 'plot_Ubreve');
%     xlabel('X'); ylabel('Y');
%     grid(axes1, 'on'); axis(axes1, 'equal'); hold(axes1, 'on');
%     % plot U, the noisy, incomplete, moving dataset
%     plot3( axes1, ...
%            U(1,:), ...
%            U(2,:), ...
%            U(3,:), ...
%            'or', ...
%            'Tag', 'plot_U');

    %%

    % prepare variable to contains all rmse
    cf  = zeros(length(r_z), length(t_z));

    loop_z = length(r_z);
    loop_t = length(t_z);

    tic;
    parfor current_z = 1:loop_z

        cf_t  = zeros(1, loop_t);

        for current_t = 1:loop_t

            % transform Ŭ with respected transformation
            U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
            scene_ptCloud = U_breve_prime';

%             % RMSE
%             [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
%             cf_t(current_t) = mean(nearest_dist);

            % GMM L2 Distance
            scale = 40e-4;
%             scale = 47.5e-4;
            [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
            cf_t(current_t) = -f;

%             % GMM CPD Distance (?)
%             X = scene_ptCloud;
%             Y = model_ptCloud;
%             sigma = 9e-8;
%             w     = 1e-90;
%             [ ~, ~, ~, negativeLogLikelihood ] = computeEStep(X, Y, sigma, w);
%             cf_t(current_t) = negativeLogLikelihood;


        end

        cf(current_z, : )  = cf_t;

    end
    toc;

    % display
    [X,Y] = meshgrid(r_z, t_z);

    figure2 = figure(2);
    surf(X,Y, cf);
    xlabel('Rz (deg)');
    ylabel('tz (mm)');
    zlabel('GMM L2 distance');
    view(90, 90);

    % search min
    minValue = min(cf(:));
    [costfunctions_min(trial, 1), costfunctions_min(trial, 2)] = find(cf == minValue);

end

middle = ceil(loop_z/2);
fprintf('rz_idx: %d\t\t\trt_idx: %d\n', costfunctions_min(1,1)-middle, costfunctions_min(1,2)-middle);
fprintf('rz_est: %.2f deg\trt_est: %.2f mm\n', r_z(costfunctions_min(:, 1)), t_z(costfunctions_min(:, 2))*1000 );

% save('results\gmm1vsgmm2.mat', 'costfunctions_min', 'r_z', 't_z');