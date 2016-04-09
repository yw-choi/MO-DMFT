        subroutine fft_1d(N,in_fft,out_fft)

        implicit none
         
        include "fftw3.f" 
   
        integer*8:: N, plan
        double precision:: in_fft(N)
        double complex:: out_fft(N/2+1)

        call dfftw_plan_dft_r2c_1d(plan,N,in_fft,out_fft,FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan,in_fft,out_fft)   ! cautious! dft_r2c! not dft
        call dfftw_destroy_plan(plan)
    
        return
        end 

        subroutine cfft_1d(N,in_fft,out_fft)
 
        implicit none
     
        include "fftw3.f"
 
        integer*8::N, plan
        double complex:: in_fft(N), out_fft(N)

        call dfftw_plan_dft_1d(plan,N,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,in_fft,out_fft)
        call dfftw_destroy_plan(plan)

        return
        end
 
        subroutine cfft_1d_bk(N,in_fft,out_fft)
 
        implicit none
     
        include "fftw3.f"
 
        integer*8::N, plan
        double complex:: in_fft(N), out_fft(N)

        call dfftw_plan_dft_1d(plan,N,in_fft,out_fft,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,in_fft,out_fft)
        call dfftw_destroy_plan(plan)

        return
        end
