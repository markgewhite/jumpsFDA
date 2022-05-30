# JumpsFDA
### Functional data analysis applied to vertical ground reaction forces obtained from countermovement jumps.

This code supports my paper entitled 'The Effects of Curve Registration on Vertical Ground Reaction Force curves of the Countermovement Jump', which is soon to be published by the Journal of Biomechanics. Link forthcoming.

The code generates 1024 linear models that predict jump height, peak power or classify the jumps into those performed with or without arm swing. Each model is based on different data preprocessing operations for:
- curve length standardisation (padding or time normalisation)
- landmark registration (16 combinations of four landmarks, including no registration)
- continuous registration (either applied or not)
- feature extraction (Functional PCA or Analysis of Characteristiing Phases)
- component rotation (unrotated/varimax)

## Key operations

The analysis is run by [code/jumpsFDA.m](code/jumpsFDA.m). The main loop iterates through the following operations:
- data partitioning
- length standardisation
- functional smoothing
- landmark registration
- continuous registration
- decomposition analysis
- linear modelling 


It depends on the Matlab library for  functional data analysis which can be found [here](https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/).


## Raw data

[jumpsFDA](code/jumpsFDA.m) reads the raw datafile [data/compactjumpsdata.mat](/data/compactjumpsdata.mat) which contains anonymised data. 
- <code>bwall: 64x16 double</code> (bodyweights)
- <code>grf.raw: 64x16 cell</code> (raw 1D VGRF data of variable lengths)
- <code>grf.initiation: 64x16 double</code> (jump initiaion index)
- <code>grf.takeoff: 64x16 double</code> (take-off index)
- <code>jumpOrder: 104x16 string</code> (master list specifying the the jump type for each trial)
- <code>nJumpsPerSubject: 64x1 double</code> (how many jumps for each subject)
- <code>sDataID: 1x64 double</code> (lookup table linking array index to subject ID)
- <code>sJumpID: 1x104 double</code> (lookup table linking trial index to jump order ID)



## Jump initiation

The countermovement jumps were performed by 55 participants, as described in the paper, and the vertical ground reaction forces recorded. These recordings started before the jump and finished after the landing. The investigation was only interested in the forces up to take-off. It was necessary to identify when the jump was started. Identifying the start and end points was the purpose of the function [demarcateJump](code/demarcateJump.m).

Take-off is easily found where VGRF drops below 10 N.

Jump initiation requires the following steps:
1.   Smooth the data with a moving average using a Gaussian window, centring it first and normalizing to bodyweight.
```Matlab
vgrf = smoothdata( (vgrf-bw)/bw, 'Gaussian', 21 );
```
2.   Detect the first significant movement when <code>threshold1</code> is breached
```Matlab
t3 = find( abs(vgrf) > threshold1, 1 );
```
3.   Work backwards to find where the VGRF falls below a smaller <code>threshold2</code>. 
```Matlab
t1 = t3;
while t1>1  && t3>1
    % find the next point backwards where VGRF is within threshold
    t1 = find( abs(vgrf(1:t3)) < threshold2, 1, 'Last' );
    if isempty(t1)
        t1 = 1;
    end

    ...

    % from there, find where the VGRF rises back above this threshold
    t3 = find( abs(vgrf(1:t1)) > threshold2, 1, 'Last' );
    if isempty(t3)
        t3 = 1;
    end

end
```
4.   Finally, find where the VGRF passes through bodyweight, which is 0 on the centred scale.
```Matlab
    % find the last zero crossover point, nearest take-off
    crossover = vgrf(1:t1).*vgrf(2:t1+1);
    t0 = find( crossover < 0, 1, 'last' );
    if isempty(t0)
        % if it does not pass through zero, find the point nearest zero
        [~, t0 ] = min( abs(vgrf(1:t1)) );
    end
    t1 = max( t0 + t3 - 1, 1 );
```

The routine also provides helpful VGRF plots if the <code>plotVGRF</code> flag is set.

The demarcation function is called by [extractVGRFData](code/extractVGRFData.m), which also excludes outliers in terms of time series length. Lenght distribution information is provided. 



## Registration

The code carries out an iterative procedure in [registerVGRF](code/registerVGRF.m), performing either landmark or continuous registration on each curve set. The purpose of registration is to more closely align the curves with one another by warping the time domain. Alignment in landmark registration is when the defined landmarks (salient features) of each curve line up with the corresponding mean landmark positions. Alignment in continuous registration attempts to align the whole length of the curve according to a specified criterion.

The following steps are performed in each registration iteration:
1.   Perform either landmark or continuous registration as specified in the setup, which returns the registered curves <code>XFdReg</code> and the warping curves <code>warpFd</code>.
2.   Update the cumulative warp in sub-function <code>updateTotalWarp</code> to keep a track of the total warping across all iterations, updating the corresponding warping function - see below.
3.   Check the cumulative warp functions are still monotonic - ie time should always go forward.
4.   Check all landmarks are properly ordered and not in extreme positions (sub-function <code>validateLandmarks</code>).
5.   Check for registration convergence indicated by the independence constant *C* (<code>decomp.c</code>) changing by less than 0.001 from the previous iteration.
6.   Store the final registered and warping curves *unless* some are invalid *and* if all must be valid, as specified in the setup.


### Landmarks

The landmarks are located by [findGRFlandmarks](code/findGRFlandmarks.m), specifically:
-  VGRF minimum
-  Peak negative power
-  Zero power (start of the concentric phase)
-  Peak positive power

### Cumulative time warping

In an iterative registration procedure it is necessary to keep a track of the total warp so the subsequent models can use the full information this provides. The sub-function <code>updateTotalWarp</code> does this by "warping the warp". The idea is to project the new time warp onto the old time warp. It works on a high-resoluition time series (1 ms intervals) rather than the functions. 

```Matlab
function totWarpT = updateTotalWarp( totWarpT, newWarpFd, t )

    N = size( totWarpT, 2 );
    for i = 1:N
        % warp the previous warp
        totWarpT(:,i) = eval_fd( totWarpT(:,i), newWarpFd(i) );
        % separate the points evenly
        totWarpT(:,i) = interp1( t, totWarpT(:,i), t, 'spline', 'extrap' );
        % impose limits in case of over/underflow
        totWarpT(:,i) = max( min( totWarpT(:,i), t(end) ), t(1) );
    end

end
```

The first step is to evaluate the newly generated warping function at the previous time series points. 

```Matlab
        totWarpT(:,i) = eval_fd( totWarpT(:,i), newWarpFd(i) );
```

If this is the first iteration the points are evenly spaced, as the function would be usually rendered. But if the points were a previous warp the warp function will be evaluated unevenly with some points close together and others further apart. 

The next step is to separate those points evenly using interpolation using the evenly spaced original linear series.

```Matlab
        totWarpT(:,i) = interp1( t, totWarpT(:,i), t, 'spline', 'extrap' );
```

This can sometimes go slightly awry at the ends of the time series where the end points can go beyond the specified limits of the time domain. The last step is to re-impose those limits.

```Matlab
        totWarpT(:,i) = max( min( totWarpT(:,i), t(end) ), t(1) );
```

### Re-smoothing the total warp

The total warp is computed as a time series, as above, but it needs to be converted back into a smooth function using the <code>resmoothWarp</code> sub-function. The re-smoothed warping functions need to be of a higher order because they may represent a more complex shapes, the compounding effect of multiple warpings. A 5th order b-spline basis function works well with a minimal roughness penalty (1E-10) to ensure high fidelity with the total warp time series, <code>warpT</code>.

```Matlab
        wBasis = create_bspline_basis( [t(1),t(end)], b, ...
                                       setup.basisOrder ); 
        wFdPar = fdPar( wBasis, ...
                        setup.basisOrder-2, ...
                        setup.wLambda );
        warpFd = smooth_basis( t, warpT, wFdPar );
```

The variable, <code>b</code>, is the number of basis functions which are adaptively increased until all function are monotonic. The loop double <code>b</code> each time or until an upper limit is reached (256 + basis order). This adaptive approach is necessary rather than just specifiying a large number of basis functions in order to save on memory/disk storage. If the maximum were used everytime, the memory required for the warping functions would be larger than all other variables. Here is the full function (note that the <code>setup</code> is a re-warping specific setup):

```Matlab
function [ warpFd, isMonotonic ] = resmoothWarp( t, warpT, setup )

    % re-smooth the warping curve using a more extensive basis
    % to ensure the closest adherence to the interpolated points
    % using the minimal number of basis functions required 

    b = 4 + setup.basisOrder;
    allMonotonic = false;
    while ~allMonotonic && b < setup.nBasisMax

        % double number of basis functions
        b = (b-setup.basisOrder)*2 + setup.basisOrder;

        % smooth time series with the trial basis
        wBasis = create_bspline_basis( [t(1),t(end)], b, ...
                                       setup.basisOrder ); 
        wFdPar = fdPar( wBasis, ...
                        setup.basisOrder-2, ...
                        setup.wLambda );
        warpFd = smooth_basis( t, warpT, wFdPar );
    
        % compute the first derivative
        warpDT = eval_fd( t, warpFd, 1 );
        % check monotonicty by curve excluding the last 5 ms
        isMonotonic = all( warpDT( 1:end-5,: ) > 0 );
        allMonotonic = all( isMonotonic );
        
    end
    disp(['Resmooth warp: number of basis functions = ' num2str(b)]);

end

```

Note that the final 5 ms before take-off are excluded for the test of monotonicity because with a very few registrations (3/128) the warping gradient turns negative at the end of the time series. Only the large dataset with CMJ<sub>NA</sub> and CMJ<sub>A</sub> are affected when continuous registration is applied to curves previously registered with 3-4 landmarks: PAD-01111C (87 curves), PAD-1111C (19) and LTN-01111C (38). The warping gradients in these curves only turn negative in the last few milliseconds. This may be due to the rapid changes in VGRF combined with the boundary effect from the end of the time series. These curves need not be excluded if the monotonicity test ignores the last few points.

 
### Validity check: VGRF curves 

Registration can very occasionally distort curves excessively. A good way to check if the VGRF curve is still reasonable is to use all landmarks irrespective of landmark registration is being performed. The sub-function <code>validateLandmarks( tSpan, XFdReg )</code> performs two checks:
- Whether a given landmark for an individual curve is out of kilter with the others in the curve set
- Whether the landmarks for each curve at in the correct order

The first check involves comparing a landmark's position with the distribution that landmark's positions in other curves. Outliers are to be expected and might be quite reasonable. A faulty registration however will typically result in the landmark's position being wildly off. Hence, extreme outliers need to be identified. Usually, the way to do is to set a Z-score threshold beyond which the outlier is considered extreme. The difficulty is that with the presence of such outliers would greatly exaggerate the SD. Therefore, the code excludes outliers outside the 95% interval (2.5%-97.5% percentile range) when calculating mean and SD.

```Matlab
% check for any extreme outliers
for i = 1:nLMs
    % compute mean and SD with outliers excluded
    limits = prctile( lm.case(:,i), [2.5 97.5] );
    outliers = (lm.case(:,i)<limits(1) | lm.case(:,i)>limits(2));
    sd = std( lm.case(~outliers, i) );
    avg = mean( lm.case(~outliers, i) );
    % check if within max Z score 
    % (very large because SD computed without outliers)
    isValid = isValid & (abs(lm.case(:,i)-avg)/sd < 50);
end
```

The second check is reasonably straightforward, checking the index positions of the sorted landmarks for a given curve.

```Matlab
% check if any landmarks are out of sequence
for i = 1:nCurves
    [~, order] = sort( lm.case(i,:) );
    isValid(i) = isValid(i) && all(diff(order)==1);
end
```


## Outputs

[jumpsFDA](code/jumpsFDA.m) outputs a Matlab data file [results/jumpsAnalysis.mat](results/jumpsAnalysis.mat). It contains the following:
- <code>decomp (2x2x16x2 cell array)</code> amplitude/phase decomposition  
- <code>fdPar (2x2 cell array)</code> FDA smoothing parameters, dataset by length standardisation
- <code>isValid (2x2x16 cell array)</code> logical array indicating valid registrations
- <code>models (2x2x16x2 cell array)</code> structures holding the fitted models
- <code>name (2x2x16x2 string array)</code> names of the datasets
- <code>vgrfACP (2x2x16x2 cell array)</code> structures holding the ACP unrotated and varimax scores
- <code>vgrfFd (2x2x16x2 cell array)</code> functional data objects for the smoothed VGRF curves
- <code>vgrfPCA (2x2x16x2 cell array)</code> structures holding the functional PCA, unrotated and varimax
- <code>warpFd (2x2x16x2 cell array)</code> functional data objects for the time-warp curves


These 4D arrays have this structure:
- 2 datasets: CMJ<sub>NA</sub> & CMJ<sub>A</sub> 
- 2 length standardisations: PAD & LTN
- 16 landmark registration: 0000-1111
- 2 continuous registrations: - & C


### Models

As 4D arrays they are not easily viewed but running the following commands compiles the results into long tables:

```Matlab
longPerformance = compileResults( models, 'perf' );
longTStat = compileResults( models, 'tStat' );
longCoeffRSq = compileResults( models, 'coeffRSq' );
longDecomp = compileResults( decomp );
faultyCountWOA = squeeze(cellfun( @sum, isValid(1,:,:) ))';
faultyCountALL = squeeze(cellfun( @sum, isValid(2,:,:) ))';
````

<code>longPerformance</code> is the most important table. It summaries the models' performance and is used to obtain the rankings presented in the paper. It has been placed in [results/ModelOutputs.xlsx](results/ModelOutputs.xlsx) for convenient analysis.


### Meta-models

The meta models can be obtained using the following command once the results file has been loaded.

```Matlab
[metaJHtov, metaJHwd, metaPP, metaCL, metaAll ] = fitMetaModels( models, 'interactions', true )
```

The second argument requests models with interactions between the predictors, as used in the paper. Alternatively, <code>'linear'</code> fits  a more straightforward model without interactions. The third argument, <code>true</code>, indicates that the predictor type should be split up into PCA/ACP and U/V (unrotated/varimax), as in the paper. The outputted variables from the function are model objects. The coefficient tables have been placed in [results/MetaModels.xlsx](results/MetaModels.xlsx) for convenience.


### Figures

Other files are dedicated to generating plots of the results. The paper focused on jump height (take-off velocity) and classification. However, models were also produced for jump height (work done) and peak power.


```Matlab
load( 'results/jumpsAnalysis.mat' )
	
[fig1, fig2 ] = plotCurves( vgrfFd, warpFd )
fig3 = plotDecompScatter( decomp )
fig4 = plotKeyComponents( vgrfPCA )
	
fig5 = plotCoeffRSq( models, ["JHtov" "jumpType"] );
```

Or, alternatively:
```Matlab
fig5 = plotCoeffRSq( models, ["JHwd" "PP"] );
```


