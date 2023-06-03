# System Requirements
## Hardware requirements
CS-FEP requires a standard computer with GPUs like NVIDIA GeForce RTX 2080 Ti or Tesla V100-SXM3-32GB.
## Software requirements
### OS Requirements
This package is supported for *Linux*. The package has been tested on the following systems:
+ Linux: CentOs 7

### Simulation Software
+ Amber with version >= 18 with self-modifying patch of the ti.F90 in AMBERHOME/src/pmemd/src in order to avoid the reorder the atoms of ligand.
- To modify: The "ti_check_and_adjust_crd" subroutine of ti.F90 
```fortran
  !reorder atoms so that if linear atoms are out of order in the prmtop
  !the simulation will still run
  !by zli: no need to do this if using our method. corresponding codes are commented.
  !ti_latm_lst_sort(:) = 0
  !do i = 1, ti_latm_cnt(1)
  !  do j = 1, ti_latm_cnt(1)
  !    atm_i = ti_latm_lst(1,i)
  !    atm_j = ti_latm_lst(2,j)
  !    total_idx = 0
  !    do m = 1, 3
  !      crd_diff = abs(crd(m, atm_i) - crd(m, atm_j))
  !      if (crd_diff .lt. 0.3d0) then
  !        total_idx = total_idx + 1
  !      end if
  !    end do
  !    if (total_idx .eq. 3) then
  !      ti_latm_lst_sort(i) = atm_j
  !    end if
  !  end do
  !end do

  !ti_latm_lst(2,:) = ti_latm_lst_sort(:)

  !do i = 1, ti_latm_cnt(1)
  !  atm_i = ti_latm_lst(1,i)
  !  atm_j = ti_latm_lst(2,i)
  !  if (atm_j .eq. 0) then
  !    write (mdout,'(a,i7,a)') '     Error : Atom ', &
  !           atm_i,' does not have match in V1 !'
  !    call mexit(mdout, 1)
  !  end if
  !end do

  !do i = 1, ti_latm_cnt(1)
  !   atm_i = ti_latm_lst(1,i)
  !   atm_j = ti_latm_lst(2,i)
  !   do m = 1, 3
  !      crd_diff = abs(crd(m, atm_i) - crd(m, atm_j))
  !      if (crd_diff .gt. 0.d0) then
  !         if (crd_diff .gt. 0.3d0) then
  !            write (mdout,'(a,i7,a,i7,a)') '     WARNING: Local coordinate ', &
  !                 atm_i,' differs from partner coordinate ', atm_j,' !'
  !            write (mdout,'(a)') &
  !               '     Atom coordinate disagreement, check input files.'
  !            call mexit(mdout, 1)
  !         else
  !            nadj = nadj + 1
  !            crd(m, atm_j) = crd(m, atm_i)
  !            if (nadj .lt. 11) then
  !               if (nadj .lt. 10) then
  !                  write (mdout,'(a,i7,a,i7,a)') &
  !                     '     WARNING: Local coordinate ', &
  !                     atm_i,' differs from partner coordinate ', atm_j,' !'
  !                  write (mdout,'(a)') &
  !                     '     Deviation is small, changing partner coordinate.'
  !               else
  !                  write (mdout,'(a)') &
  !                     '     ... making more adjustments ...'
  !               end if
  !            end if
  !         end if
  !      end if
  !   end do
  !end do
```
The installation guide of Amber: http://ambermd.org/Installation.php
- Typical install time of Amber on a "normal" desktop computer: about 40 mins.
### Analysis Code Dependencies
- The analysis code depends on the Python scientific stack. 
```
numpy
pandas
pymbar
matplotlib
```
- Typical install time of python and required packages on a "normal" desktop computer: about 20 mins.


# Running simulation
- Use this script to run the MD simulation
`run.sh`
```sh
#!/bin/bash
ROOT=`pwd`
TARGET='cdk2_28-26_demo'
WINS='0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00'
SIDES='complex/cM2A complex/cM2B ligands/lM2A ligands/lM2B'
_CUDA_VISIBLE_DEVICES_=0
for side in $SIDES;do
    cd $ROOT/$TARGET/$side
    for win in $WINS;do
        cd $win
        sh submit.sh $_CUDA_VISIBLE_DEVICES_
    done
done
```
- Expected run time for demo's MD simulation on a "normal" desktop computer(with one gpu like Tesla V100-SXM3-32GB): 6-8 hours


# Analysis
- The demo's simulation out was provided in the cdk2_28-26.tar.bz2, use the following command to extract the data:
```sh
tar -zxvf cdk2_28-26.tar.bz2
```
- Use the in-house free energy calculation python script: `analysis_free_energy.py`
- `aly.sh`
```sh
#!/bin/bash
ROOT=`pwd`
FEP_PATH="$ROOT/cdk2_28-26/FEP"
mkdir -p analysis_out
cd analysis_out
python $ROOT/analysis_free_energy.py -d ${FEP_PATH} -f 0.25 -o fe_out.csv
cd ..
```
- Expected run time for demo's analysis on a "normal" desktop computer: 3-5 mins
