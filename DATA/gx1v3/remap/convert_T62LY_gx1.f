!***********************************************************************

      program convert_T62_gx1

!coyote:
!ifort convert_T62LY_gx1.f -o convert_gx1 -w95 -r8 -i4 -assume byterecl -convert big_endian -I/usr/include -I/usr/projects/climate/maltrud/local/include_coyote -L/usr/projects/climate/maltrud/local/lib_coyote -lnetcdf
!-----------------------------------------------------------------------
!
!     this program remaps Large/Yeager atmosphere data from a T62 grid to 
!     the gx1 ocean grid
!
!-----------------------------------------------------------------------

      implicit none

      integer :: yr

      do yr = 1980,1980  ! date in filenames below:  2004_08_30
!      do yr = 1981,2000  ! date in filenames below:  2004_08_03
!      do yr = 2001,2004  ! date in filenames below:  2006_04_06
!      do yr = 2005,2006  ! date in filenames below:  7JULY2008
        call grid_conversion (yr)
      enddo

      end program convert_T62_gx1

!***********************************************************************

      subroutine grid_conversion (year)

      implicit none
      
      include 'netcdf.inc'

!-----------------------------------------------------------------------
!
!     variables to be mapped
!
!-----------------------------------------------------------------------

      integer year, n, out_rec

      integer, dimension(3) ::
     &  beg,       ! beginning point for netCDF array
     &  cnt        ! count for each axis of netCDF array 

      integer ::
     &  ncid_u_10,
     &  ncid_v_10,
     &  ncid_t_10,
     &  ncid_q_10

      integer, parameter ::
     &  unit_u_10  = 21,
     &  unit_v_10  = 22,
     &  unit_t_10  = 23,
     &  unit_q_10  = 24

      character(60) :: atm_file 

      character(60), parameter :: 
     &  remap_file = 'rmp_T62_to_gx1v3_bil.nc'

      character(60) :: 
     &  out_file_t_10,
     &  out_file_q_10,
     &  out_file_u_10,
     &  out_file_v_10

      real, dimension(:), allocatable ::
     &      atm_field

!      double precision, dimension(:), allocatable ::
      real, dimension(:), allocatable ::
     &      ocn_field  ! remapped field to be written

!-----------------------------------------------------------------------
!
!     variables required for mapping
!
!-----------------------------------------------------------------------

      integer :: ! netCDF ids
     &         ncstat, nc_file_id, nc_time_id,
     &         nc_srcgrdsize_id, nc_dstgrdsize_id,
     &         nc_srcgrddims_id, nc_dstgrddims_id,
     &         nc_srcgrdrank_id, nc_dstgrdrank_id,
     &         nc_srcfield_id, nc_dstfield_id,
     &         nc_numlinks_id, nc_numwgts_id, 
     &         nc_srcgrdadd_id, nc_dstgrdadd_id, nc_rmpmatrix_id

      integer, dimension(3) ::
     &   nc_srcdims3_id,  ! netCDF ids for 2d array dims
     &   nc_dstdims3_id,  ! netCDF ids for 2d array dims
     &   istart,          ! starting address for array block
     &   icounta,         ! count along each dimension for block
     &   icounto          ! count along each dimension for block

      integer link, num_wts, num_links, ocn_size, atm_size, 
     &                                  ocn_rank, atm_rank

      integer, dimension(:), allocatable ::
     &         ocn_add, atm_add,     ! addresses for remapping
     &         ocn_dim, atm_dim      ! dimensions of each axis

      double precision, dimension(:,:), allocatable ::
     &         wgts

!-----------------------------------------------------------------------
!
!     open remap file
!
!-----------------------------------------------------------------------

      print *,'opening remap file'
      ncstat = nf_open(remap_file, NF_NOWRITE, nc_file_id)

!-----------------------------------------------------------------------
!
!     read dimension information
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_size', 
     &                      nc_srcgrdsize_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdsize_id, atm_size)
      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_size', 
     &                      nc_dstgrdsize_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdsize_id, ocn_size)
      ncstat = nf_inq_dimid(nc_file_id, 'num_links', 
     &                      nc_numlinks_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, 
     &                       num_links)
      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_rank', 
     &                      nc_srcgrdrank_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdrank_id, atm_rank)
      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_rank', 
     &                      nc_dstgrdrank_id)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdrank_id, ocn_rank)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      print *,' ocn_size = ',ocn_size,' atm_size = ',atm_size

      allocate( ocn_field(ocn_size),
     &          atm_field(atm_size), 
     &          ocn_dim  (ocn_rank),
     &          atm_dim  (atm_rank),
     &          ocn_add  (num_links),
     &          atm_add  (num_links),
     &          wgts     (num_wts,num_links))

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'src_address', 
     &                                   nc_srcgrdadd_id)
      ncstat = nf_inq_varid(nc_file_id, 'dst_address', 
     &                                   nc_dstgrdadd_id)
      ncstat = nf_inq_varid(nc_file_id, 'remap_matrix', 
     &                                   nc_rmpmatrix_id)
      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, 
     &                        atm_add)
      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, 
     &                        ocn_add)
      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                                       wgts)
      ncstat = nf_inq_varid(nc_file_id, 'src_grid_dims', 
     &                                   nc_srcgrddims_id)
      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_dims', 
     &                                   nc_dstgrddims_id)
      ncstat = nf_get_var_int(nc_file_id, nc_srcgrddims_id, 
     &                        atm_dim)
      ncstat = nf_get_var_int(nc_file_id, nc_dstgrddims_id, 
     &                        ocn_dim)

!-----------------------------------------------------------------------
!
!     close remap file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)

      print *,'finished reading remap data'

!-----------------------------------------------------------------------
!
!     open output files
!
!-----------------------------------------------------------------------


      write (out_file_t_10,1000), 't_10.',year
      write (out_file_q_10,1000), 'q_10.',year
      write (out_file_u_10,1000), 'u_10.',year
      write (out_file_v_10,1000), 'v_10.',year
 1000   format('DATA/gx1v3/LargeYeager/4XDAILY/',a5,i4,'.dat')

      open(unit_u_10,     file=out_file_u_10, status='unknown',
     &             form='unformatted', access='direct', recl=8*ocn_size)
      open(unit_v_10,     file=out_file_v_10, status='unknown',
     &             form='unformatted', access='direct', recl=8*ocn_size)
      open(unit_t_10,     file=out_file_t_10, status='unknown',
     &             form='unformatted', access='direct', recl=8*ocn_size)
      open(unit_q_10,  file=out_file_q_10, status='unknown',
     &             form='unformatted', access='direct', recl=8*ocn_size)

!-----------------------------------------------------------------------
!
!     read in data and remap
!
!-----------------------------------------------------------------------

      beg(1) = 1
      beg(2) = 1
      cnt(1) = atm_dim(1)
      cnt(2) = atm_dim(2)
      cnt(3) = 1

        atm_file = ' '
        write(atm_file,1001) year
 1001   format('T62_LY/u_10.',i4,'.2004_08_30.nc')
! 1001   format('T62_LY/u_10.',i4,'.2006_04_06.nc')
! 1001   format('T62_LY/u_10.',i4,'.7JULY2008.nc')
        print *,'reading file ',atm_file
        ncstat = nf_open(atm_file,NF_NOWRITE,nc_file_id)
        ncstat = nf_inq_varid(nc_file_id,'U_10_MOD' ,ncid_u_10)
        out_rec = 0
        do n=1,1460
!          print *,'u_10',n
          out_rec = out_rec + 1
          beg(3) = n
          ncstat = nf_get_vara_double(nc_file_id,ncid_u_10,beg,cnt,
     &                                atm_field)
          ocn_field = 0.0
          do link=1,num_links
            ocn_field(ocn_add(link)) = ocn_field(ocn_add(link)) +
     &                    wgts(1,link)*atm_field(atm_add(link))
          end do
          write(unit_u_10,rec=out_rec) ocn_field

          if (n.eq.1) then
          print *,minval(atm_field),maxval(atm_field),
     &                    sum(atm_field)/atm_size
          print *,minval(ocn_field),maxval(ocn_field),
     &                    sum(ocn_field)/ocn_size
          endif
        end do
        ncstat = nf_close(nc_file_id)

        atm_file = ' '
        write(atm_file,1002) year
! 1002   format('T62_LY/v_10.',i4,'.2004_08_03.nc')
! 1002   format('T62_LY/v_10.',i4,'.2006_04_06.nc')
 1002   format('T62_LY/v_10.',i4,'.7JULY2008.nc')
        print *,'reading file ',atm_file
        ncstat = nf_open(atm_file,NF_NOWRITE,nc_file_id)
        ncstat = nf_inq_varid(nc_file_id,'V_10_MOD' ,ncid_v_10)
        out_rec = 0
        do n=1,1460
!          print *,'v_10',n
          out_rec = out_rec + 1
          beg(3) = n
          ncstat = nf_get_vara_double(nc_file_id,ncid_v_10,beg,cnt,
     &                                atm_field)
          ocn_field = 0.0
          do link=1,num_links
            ocn_field(ocn_add(link)) = ocn_field(ocn_add(link)) +
     &                    wgts(1,link)*atm_field(atm_add(link))
          end do
          write(unit_v_10,rec=out_rec) ocn_field

          if (n.eq.1) then
          print *,minval(atm_field),maxval(atm_field),
     &                    sum(atm_field)/atm_size
          print *,minval(ocn_field),maxval(ocn_field),
     &                    sum(ocn_field)/ocn_size
          endif
        end do
        ncstat = nf_close(nc_file_id)

        atm_file = ' '
        write(atm_file,1003) year
! 1003   format('T62_LY/t_10.',i4,'.2004_08_03.nc')
! 1003   format('T62_LY/t_10.',i4,'.2006_04_06.nc')
 1003   format('T62_LY/t_10.',i4,'.7JULY2008.nc')
        print *,'reading file ',atm_file
        ncstat = nf_open(atm_file,NF_NOWRITE,nc_file_id)
        ncstat = nf_inq_varid(nc_file_id,'T_10_MOD' ,ncid_t_10)
        out_rec = 0
        do n=1,1460
!          print *,'t_10',n
          out_rec = out_rec + 1
          beg(3) = n
          ncstat = nf_get_vara_double(nc_file_id,ncid_t_10,beg,cnt,
     &                                atm_field)
          ocn_field = 0.0
          do link=1,num_links
            ocn_field(ocn_add(link)) = ocn_field(ocn_add(link)) +
     &                    wgts(1,link)*atm_field(atm_add(link))
          end do
          write(unit_t_10,rec=out_rec) ocn_field

          if (n.eq.1) then
          print *,minval(atm_field),maxval(atm_field),
     &                    sum(atm_field)/atm_size
          print *,minval(ocn_field),maxval(ocn_field),
     &                    sum(ocn_field)/ocn_size
          endif
        end do
        ncstat = nf_close(nc_file_id)

        atm_file = ' '
        write(atm_file,1004) year
! 1004   format('T62_LY/q_10.',i4,'.2004_08_03.nc')
! 1004   format('T62_LY/q_10.',i4,'.2006_04_06.nc')
 1004   format('T62_LY/q_10.',i4,'.7JULY2008.nc')
        print *,'reading file ',atm_file
        ncstat = nf_open(atm_file,NF_NOWRITE,nc_file_id)
        ncstat = nf_inq_varid(nc_file_id,'Q_10_MOD' ,ncid_q_10)
        out_rec = 0
        do n=1,1460
!          print *,'q_10',n
          out_rec = out_rec + 1
          beg(3) = n
          ncstat = nf_get_vara_double(nc_file_id,ncid_q_10,beg,cnt,
     &                                atm_field)
          ocn_field = 0.0
          do link=1,num_links
            ocn_field(ocn_add(link)) = ocn_field(ocn_add(link)) +
     &                    wgts(1,link)*atm_field(atm_add(link))
          end do
          write(unit_q_10,rec=out_rec) ocn_field

          if (n.eq.1) then
          print *,minval(atm_field),maxval(atm_field),
     &                    sum(atm_field)/atm_size
          print *,minval(ocn_field),maxval(ocn_field),
     &                    sum(ocn_field)/ocn_size
          endif
        end do
        ncstat = nf_close(nc_file_id)

!-----------------------------------------------------------------------
!
!     close input file
!
!-----------------------------------------------------------------------

      deallocate( ocn_field, atm_field, ocn_add, atm_add, wgts)
      close(unit_u_10)
      close(unit_v_10)
      close(unit_t_10)
      close(unit_q_10)

!-----------------------------------------------------------------------

      end subroutine grid_conversion

!***********************************************************************

