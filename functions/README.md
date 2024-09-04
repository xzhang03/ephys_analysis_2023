## Function list

#### Load abf files (not written by me)
```Matlab
[d,si,h]=abfload(fn,varargin)
```

#### Calculate charge transfer from V clamps
```Matlab
ephysChargeTransfer(fpath, varargin)
```

#### Calculate current sizes from EPSC/IPSC events
```Matlab
outputmat = ephysCurrent(fpath, varargin)
```

#### A UI to find EPSCs
```Matlab
fn_out = ephysEPSCFinder(fpath, varargin)
```

#### A UI to find IPSCs
```Matlab
fn_out = ephysIPSCFinder(fpath, varargin)
```

#### Calculate instantaneous firing rates from loose patch data
```Matlab
outputmat = ephysInstaFR(fpath, varargin)
```

#### A UI to find spikes from loose-patch data
```Matlab
fn_out = ephysLoosePatchSpikeFinder(fpath, varargin)
```

#### Find when stims are done from ephys digitized data
```Matlab
fn_out = ephysStimFinder(fpath, varargin)
```

#### Calculate moving percentiles
```Matlab
outvec = movprctile(input, window, percentile)
```


#### 
```Matlab

```
