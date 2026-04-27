C=======================================================================
C  generate_dis.f
C  DIS-only generator for SIMC: read events from ASCII dis_datfile
C
C  Expected columns (SIMC units):
C    EvtID   pe_px  pe_py  pe_pz   pr_px  pr_py  pr_pz   pidrec   zvtx   weight
C
C  Units:
C    momentum components: MeV/c
C    zvtx: cm   (relative to target center; targ%zoffset will be added)
C    weight: arbitrary (user-defined); used as main%gen_weight
C
C  Output:
C    fills main%target, vertex (true), orig (=vertex), success
C=======================================================================

      subroutine generate_dis(main, vertex, orig, success)

      USE structureModule
      implicit none
      include 'simulate.inc'
C      include 'constants.inc'

      type(event_main) :: main
      type(event)      :: vertex, orig
      logical          :: success

C----------------------------
C  Persistent I/O state
C----------------------------
      integer, save :: iu = 97
      logical, save :: opened = .false.
      integer :: ios
      character*512 :: line

C----------------------------
C  Input fields
C----------------------------
      integer*8 evtid
      real*8 pe_px, pe_py, pe_pz
      real*8 pr_px, pr_py, pr_pz
      integer pidrec
      real*8 zvtx, w_in

C----------------------------
C  Local kinematics
C----------------------------
      real*8 peP, prP
      real*8 Ein, Ee, Ep
      real*8 uex, uey, uez
      real*8 upx, upy, upz
      real*8 tmp

C----------------------------
C  Constants
C----------------------------
      real*8 tiny
      parameter (tiny = 1.0d-12)

C  Fixed incident beam energy for this DIS sample: 10.2 GeV = 10200 MeV
      real*8 Ein_fixed
      parameter (Ein_fixed = 10200.d0)

      success = .false.

C===========================================================
C  Open file once
C===========================================================
      if (.not. opened) then
         if (len_trim(dis_datfile) .eq. 0) then
            write(6,*) 'generate_dis: dis_datfile is empty.'
            stop
         endif

         open(unit=iu, file=trim(dis_datfile), status='old',
     >        form='formatted', iostat=ios)

         if (ios .ne. 0) then
            write(6,*) 'generate_dis: cannot open file: ', trim(dis_datfile)
            write(6,*) 'generate_dis: iostat = ', ios
            stop
         endif

         write(6,*) 'generate_dis: opened file = ', trim(dis_datfile)
         opened = .true.
      endif

C===========================================================
C  Read next valid line (skip blank/comment/bad line)
C===========================================================
   10 continue
      read(iu,'(A)',iostat=ios) line
      if (ios .ne. 0) then
         write(6,*) 'generate_dis: reached EOF or read error (iostat=',
     >              ios, ')'
         stop
      endif

      if (len_trim(line) .eq. 0) goto 10
      if (line(1:1) .eq. '#') goto 10

      read(line,*,iostat=ios) evtid,
     >     pe_px, pe_py, pe_pz,
     >     pr_px, pr_py, pr_pz,
     >     pidrec, zvtx, w_in

      if (ios .ne. 0) goto 10

      ntup%dis_event_id = dble(evtid)

C===========================================================
C  Progress report every 100 event IDs
C===========================================================
      if (mod(evtid,100_8) .eq. 0_8) then
         write(6,*) 'generate_dis: now processing event_id = ', evtid
      endif

C Optional warning
      if (pidrec .ne. 2212) then
         if (debug(1)) then
            write(6,*) 'generate_dis: warning, pidrec != 2212 : ', pidrec
         endif
      endif

C===========================================================
C  Initialize target quantities for this event
C  (use averaged values to stay consistent with SIMC MC/recon)
C===========================================================
      main%target%x = 0.d0
      main%target%y = 0.d0
      main%target%z = zvtx + targ%zoffset
      main%target%rasterx = 0.d0
      main%target%rastery = 0.d0

      main%target%Eloss(1) = targ%Eloss(1)%ave
      main%target%Eloss(2) = targ%Eloss(2)%ave
      main%target%Eloss(3) = targ%Eloss(3)%ave

      main%target%teff(1) = targ%teff(1)%ave
      main%target%teff(2) = targ%teff(2)%ave
      main%target%teff(3) = targ%teff(3)%ave

      main%target%Coulomb = targ%Coulomb%ave

C===========================================================
C  Fixed incident beam energy: 10.2 GeV
C===========================================================
      Ein = Ein_fixed
      vertex%Ein = Ein
      main%Ein_shift = 0.d0
      main%Ee_shift  = 0.d0

C===========================================================
C  Scattered electron: from momentum components
C===========================================================
      peP = sqrt(pe_px*pe_px + pe_py*pe_py + pe_pz*pe_pz)
      if (peP .le. tiny) goto 10

      Ee  = peP
      uex = pe_px / peP
      uey = pe_py / peP
      uez = pe_pz / peP

      vertex%e%P = peP
      vertex%e%E = Ee

      tmp = max(-1.d0, min(1.d0, uez))
      vertex%e%theta = acos(tmp)
      vertex%e%phi   = atan2(uey, uex)
      if (vertex%e%phi .lt. 0.d0) vertex%e%phi = vertex%e%phi + 2.d0*pi

C Convert to spectrometer variables (xptar/yptar) consistent with SIMC
      call spectrometer_angles(spec%e%theta, spec%e%phi,
     >     vertex%e%xptar, vertex%e%yptar,
     >     vertex%e%theta, vertex%e%phi)

      vertex%e%delta = 100.d0 * (vertex%e%P - spec%e%P) / spec%e%P

C===========================================================
C  Recoil proton: assume Mp mass, momentum from components
C===========================================================
      prP = sqrt(pr_px*pr_px + pr_py*pr_py + pr_pz*pr_pz)
      if (prP .le. tiny) goto 10

      Ep  = sqrt(prP*prP + Mp*Mp)
      upx = pr_px / prP
      upy = pr_py / prP
      upz = pr_pz / prP

      vertex%p%P = prP
      vertex%p%E = Ep

      tmp = max(-1.d0, min(1.d0, upz))
      vertex%p%theta = acos(tmp)
      vertex%p%phi   = atan2(upy, upx)
      if (vertex%p%phi .lt. 0.d0) vertex%p%phi = vertex%p%phi + 2.d0*pi

      call spectrometer_angles(spec%p%theta, spec%p%phi,
     >     vertex%p%xptar, vertex%p%yptar,
     >     vertex%p%theta, vertex%p%phi)

      vertex%p%delta = 100.d0 * (vertex%p%P - spec%p%P) / spec%p%P

C===========================================================
C  Fill a few commonly used DIS kinematics
C===========================================================
      vertex%nu = vertex%Ein - vertex%e%E
      vertex%Q2 = 2.d0 * vertex%Ein * vertex%e%E *
     >            (1.d0 - cos(vertex%e%theta))
      vertex%q  = sqrt(max(0.d0, vertex%Q2 + vertex%nu*vertex%nu))

      if (vertex%nu .gt. tiny) then
         vertex%xbj = vertex%Q2 / (2.d0 * Mp * vertex%nu)
      else
         vertex%xbj = 0.d0
      endif

C Unit vector uq (q-hat), with beam along +z
      if (vertex%q .gt. tiny) then
         vertex%uq%x = -vertex%e%E * sin(vertex%e%theta) *
     >                 cos(vertex%e%phi) / vertex%q
         vertex%uq%y = -vertex%e%E * sin(vertex%e%theta) *
     >                 sin(vertex%e%phi) / vertex%q
         vertex%uq%z = (vertex%Ein - vertex%e%E *
     >                 cos(vertex%e%theta)) / vertex%q
      else
         vertex%uq%x = 0.d0
         vertex%uq%y = 0.d0
         vertex%uq%z = 1.d0
      endif

C===========================================================
C  DIS diagnostic missing quantities
C===========================================================
      vertex%Pmx = vertex%p%P * sin(vertex%p%theta) *
     >             cos(vertex%p%phi) - vertex%q * vertex%uq%x
      vertex%Pmy = vertex%p%P * sin(vertex%p%theta) *
     >             sin(vertex%p%phi) - vertex%q * vertex%uq%y
      vertex%Pmz = vertex%p%P * cos(vertex%p%theta) -
     >             vertex%q * vertex%uq%z

      vertex%Pmiss = sqrt(vertex%Pmx*vertex%Pmx +
     >                    vertex%Pmy*vertex%Pmy +
     >                    vertex%Pmz*vertex%Pmz)

      vertex%Emiss = vertex%nu + targ%M - vertex%p%E

      vertex%Pm   = vertex%Pmiss
      vertex%Em   = vertex%Emiss
      vertex%Trec = 0.d0
      vertex%Mrec = 0.d0

C===========================================================
C  Weights: keep it simple; complete_main will recompute main%weight
C===========================================================
      main%gen_weight = w_in
      main%jacobian   = 1.d0
      main%SF_weight  = 1.d0
      main%sigcc      = 1.d0

C===========================================================
C  Copy and return
C===========================================================
      orig = vertex
      success = .true.
      return

      end
