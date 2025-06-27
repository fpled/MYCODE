function [T,idSnap] = setTimeSchemeSingleEdgeCrack(Dim,symmetry,loading,test)
% function [T,idSnap] = setTimeSchemeSingleEdgeCrack(Dim,symmetry,loading,test)
%
% [T,idSnap] = setTimeSchemeSingleEdgeCrack(Dim,symmetry,loading,test)
% Returns the TIMEMODEL 'T' given the spatial dimension 'Dim', the material
% 'symmetry' the mechanical 'loading' and wether it is 'test' case or not.

switch lower(symmetry)
    case 'isotropic' % isotropic material
        if Dim==2
            switch lower(loading)
                case 'tension'
                    % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Ambati, Gerasimov, De Lorenzis, 2015, CM]
                    % du = 1e-5 mm during the first 500 time steps (up to u = 5e-3 mm)
                    % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                    dt0 = 1e-8;
                    nt0 = 500;
                    dt1 = 1e-9;
                    nt1 = 1300;
                    if test
                        dt0 = 1e-7;
                        nt0 = 50;
                        dt1 = 1e-8;
                        nt1 = 400;
                    end
                    t0 = linspace(dt0,nt0*dt0,nt0);
                    t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                    t = [t0,t1];

                    % [Miehe, Welschinger, Hofacker, 2010 IJNME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                    % du = 1e-5 mm during 630 time steps (up to u = 6.3e-3 mm)
                    % dt = 1e-8;
                    % nt = 630;
                    % t = linspace(dt,nt*dt,nt);

                    % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US]
                    % du = 1e-4 mm during 63 time steps (up to u = 6.3e-3 mm)
                    % dt = 1e-7;
                    % nt = 63;
                    % t = linspace(dt,nt*dt,nt);

                    % [Liu, Li, Msekh, Zuo, 2016, CMS]
                    % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                    % du = 1e-6 mm during the last 1300 time steps (up to u = 6.3e-3 mm)
                    % dt0 = 1e-7;
                    % nt0 = 50;
                    % dt1 = 1e-9;
                    % nt1 = 1300;
                    % t0 = linspace(dt0,nt0*dt0,nt0);
                    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                    % t = [t0,t1];
                case 'shear'
                    % [Miehe, Welschinger, Hofacker, 2010 IJNME]
                    % du = 1e-4 mm during the first 100 time steps (up to u = 10e-3 mm)
                    % du = 1e-6 mm during the last 10 000 time steps (up to u = 20e-3 mm)
                    % dt0 = 1e-7;
                    % nt0 = 100;
                    % dt1 = 1e-9;
                    % nt1 = 10000;
                    % t0 = linspace(dt0,nt0*dt0,nt0);
                    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                    % t = [t0,t1];

                    % [Liu, Li, Msekh, Zuo, 2016, CMS]
                    % du = 1e-4 mm during the first 50 time steps (up to u = 5e-3 mm)
                    % du = 1e-5 mm during the last 1500 time steps (up to u = 20e-3 mm)
                    % dt0 = 1e-7;
                    % nt0 = 50;
                    % dt1 = 1e-8;
                    % nt1 = 1500;
                    % t0 = linspace(dt0,nt0*dt0,nt0);
                    % t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                    % t = [t0,t1];

                    % [Miehe, Hofacker, Welschinger, 2010, CMAME], [Nguyen, Yvonnet, Zhu, Bornert, Chateau, 2015, EFM]
                    % [Ambati, Gerasimov, De Lorenzis, 2015, CM], [Wu, Nguyen, Nguyen, Sutula, Bordas, Sinaie, 2020, AAM],
                    % [Ulloa, Rodriguez, Samaniego, Samaniego, 2019, US], [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                    % du = 1e-5 mm during 1500 time steps (up to u = 15e-3 mm)
                    dt = 1e-8;
                    nt = 1500;
                    % nt = 2000;
                    if test
                        dt = 5e-8;
                        % nt = 300;
                        nt = 400;
                        dt = 5e-7;
                        nt = 40;
                    end
                    t = linspace(dt,nt*dt,nt);
            end
        elseif Dim==3
            dt = 1e-8;
            nt = 2500;
            if test
                dt = 1e-7;
                nt = 250;
            end
            t = linspace(dt,nt*dt,nt);
        end

    case 'anisotropic' % anisotropic material
        if Dim==2
            switch lower(loading)
                case 'tension'
                    % [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]
                    % du = 6e-5 mm during the first 200 time steps (up to u = 12e-3 mm)
                    % du = 2e-5 mm during the last 600 time steps (up to u = 24e-3 mm)
                    dt0 = 6e-8;
                    nt0 = 200;
                    dt1 = 2e-8;
                    nt1 = 600;
                    if test
                        dt0 = 6e-7;
                        nt0 = 20;
                        dt1 = 2e-7;
                        nt1 = 50;
                    end

                case 'shear'
                    % du = 1e-4 mm during the first 200 time steps (up to u = 20e-3 mm)
                    % du = 2e-5 mm during the last 2000 time steps (up to u = 60e-3 mm)
                    dt0 = 1e-7;
                    nt0 = 200;
                    dt1 = 2e-8;
                    nt1 = 2000;
                    if test
                        dt0 = 1e-6;
                        nt0 = 20;
                        dt1 = 2e-7;
                        nt1 = 200;
                    end
            end
            t0 = linspace(dt0,nt0*dt0,nt0);
            t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
            t = [t0,t1];

            if isotropicTest
                switch lower(loading)
                    case 'tension'
                        dt0 = 1e-8;
                        nt0 = 500;
                        dt1 = 1e-9;
                        nt1 = 1300;
                        if test
                            dt0 = 1e-7;
                            nt0 = 50;
                            dt1 = 1e-8;
                            nt1 = 400;
                        end
                        t0 = linspace(dt0,nt0*dt0,nt0);
                        t1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
                        t = [t0,t1];
                    case 'shear'
                        dt = 1e-8;
                        nt = 1500;
                        if test
                            dt = 5e-8;
                            nt = 400;
                        end
                        t = linspace(dt,nt*dt,nt);
                end
            end

        elseif Dim==3
            dt = 1e-8;
            nt = 2500;
            if test
                dt = 1e-7;
                nt = 250;
            end
            t = linspace(dt,nt*dt,nt);
        end
    otherwise
        error('Wrong material symmetry class');
end
T = TIMEMODEL(t);

%% Set indices of snapshots of interest to display
[t,~] = gettevol(T);
switch lower(symmetry)
    case 'isotropic'
        switch lower(loading)
            case 'tension'
                idSnap = find(abs(t-5.5e-6)<eps | abs(t-5.75e-5)<eps | abs(t-6e-6)<eps | abs(t-6.25e-6)<eps);
            case 'shear'
                idSnap = find(abs(t-1e-5)<eps | abs(t-1.25e-5)<eps | abs(t-1.35e-5)<eps | abs(t-1.5e-5)<eps);
        end
    case 'anisotropic'
        switch lower(loading)
            case 'tension'
                idSnap = find(abs(t-9e-6)<eps | abs(t-12e-6)<eps | abs(t-13.5e-6)<eps | abs(t-15e-6)<eps | abs(t-20e-6)<eps);
            case 'shear'
                idSnap = find(abs(t-20e-6)<eps | abs(t-30e-6)<eps | abs(t-40e-6)<eps | abs(t-50e-6)<eps);
        end
end
% idSnap = [idSnap,length(T)]; % indices of snapshots to display
idSnap = length(T); % index to display only last snapshot