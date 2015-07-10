#ifndef _INIT_H
#define _INIT_H

string optionsname = "options";

char texthelp[]="Usage: ./gen.x [OPTION]... [FILE]\n"
        "Generates configurations of the U(1) gauge-Higgs model with 2 flavors\n"
        "(in the dual representation) using the surface worm algorithm\n"
        "\n"
        "Mandatory arguments to long options are mandatory for short options too.\n"
        "  -l, --lambda LAMBDA    coupling of the phi^4 term (default = 1),\n"
        "  -K, --kappa KAPPA      hopping parameter,\n"
        "  -k, --dk DKAPPA        step size for kappa,\n"
        "  -B, --beta BETA        gauge coupling,\n"
        "  -b, --db DBETA         step size for beta,\n"
        "  -M, --mu MU            chemical potential (default = 0),\n"
        "  -m, --dm DMU               step size for mu (default = 0),\n"

        "  -n, --meas MEASUREMENTS   number of measurements,\n"
        "  -s, --skip SKIPS          discarded steps between meas. (default = 10),\n"
        "  -e, --equi   EQUILIBRATION    equilibration steps (default = 10^6),\n"

        "  -r, --read FILE       read configuration file FILE (default = FALSE),\n"

        "  -h  --help            display this help and exit\n"
        "\n"
        "Exit status:\n"
        " 0  if OK,\n"
        " 1  if ERROR,\n"
        "\n"
        "Report bugs to ydelgado83@gmail.com\n";

//_________________________________________________________________________
inline void read_params(int argc, char *argv[] )
{
    if ( id==0 ) cout << "Init: Initializing... " << endl;

    if(argc<2)
    {
        if ( id==0 ) cout << endl << texthelp << endl;
        return;
    }

  // Default values
  lambda = 1;
  nequi = 1000000;
  nskip = 10;
  mu = 0;
  dmu = 0;
  rconf = false;

  int c;
       
    while (1)
    {
        static struct option long_options[] =
        {
            /* These options don't set a flag.
            We distinguish them by their indices. */
            {"lambda", required_argument, 0, 'l'},
            {"kappa", required_argument, 0, 'K'},
            {"dk", required_argument, 0, 'k'},
            {"beta", required_argument, 0, 'B'},
            {"db", required_argument, 0, 'b'},
            {"mu", required_argument, 0, 'M'},
            {"dm", required_argument, 0, 'm'},
            {"meas", required_argument, 0, 'n'},
            {"skip", required_argument, 0, 's'},
            {"equi", required_argument, 0, 'e'},
            {"read", no_argument, 0, 'r'},
            /* These options set a flag. */
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
            
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "l:K:k:B:b:M:m:n:s:e:r:h",
        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");

            case 'l':
                lambda = atof(optarg);
                break;

            case 'B':
                beta0 = atof(optarg);
                break;

            case 'b':
                dbeta = atof(optarg);
                break;

            case 'K':
                kappa0 = atof(optarg);
                break;

            case 'k':
                dkappa = atof(optarg);
                break;

            case 'M':
                mu0 = atof(optarg);
                break;

            case 'm':
                dmu = atof(optarg);
                break;

            case 'n':
                nmeas = atoi(optarg);
                break;

            case 'e':
                nequi = atoi(optarg);
                break;

            case 's':
                nskip = atoi(optarg);
                break;

            case 'r':
                rconf = true;
                sprintf(conffile,"%s",optarg);
                break;

            case 'h':
                if (id == 0) cout << endl << texthelp << endl;
                abort();

            default:
                if (id == 0) cout << endl << texthelp << endl;
                abort();
        }
    }

    // print options
    if (id==0){
      cout << "lambda: " << lambda << endl;
        cout << "beta0: " << beta0 << " - dbeta: " << dbeta << endl;
        cout << "kappa0: " << kappa0 << " - dkappa: " << dkappa << endl;
        cout << "mu0: " << mu0 << " - dmu: " << dmu << endl;
        cout << "nmeas: " << nmeas << endl;
        cout << "nequi: " << nequi << endl;
        cout << "nskip: " << nskip << endl;
        if (rconf) cout << "reading confs. from: \"" << conffile << "\"" << endl;   
    }

    /* Print any remaining command line arguments (not options). */
    /* Get the name of the outputfile */
    string tmpstring;
    if ( optind >= argc ) 
    {
        if ( id==0 ) 
            cout << "Init.ERROR: file with output file..." << endl;
        return;
    }
    else tmpstring = argv[optind];

    sprintf( outfile, "%s_%d.out", tmpstring.c_str( ), id );
    cout << "outfile : \"" << outfile << "\"" << endl;

  file.open(outfile,ios::trunc | ios::out | ios::binary );

  if (id == 0 )
  {
    file.write((char*)&leng,sizeof(int));
    file.write((char*)&sub_leng,sizeof(int));
    file.write((char*)&leng_t,sizeof(int));
    file.write((char*)&nproc, sizeof(int));

    file.write((char*)&lambda,sizeof(double));

    file.write((char*)&beta0,sizeof(double));
      file.write((char*)&dbeta,sizeof(double));

    file.write((char*)&kappa0,sizeof(double));
    file.write((char*)&dkappa,sizeof(double));

    file.write((char*)&mu0,sizeof(double));
    file.write((char*)&dmu,sizeof(double));

        file.write((char*)&nequi,sizeof(int));
    file.write((char*)&nmeas,sizeof(int));
      file.write((char*)&nskip,sizeof(int));

        cout << "Init: done!\n" << endl;
  }

    return;
}

//_________________________________________________________________________
inline void read_conf( )
{
  cout << "\nreading configuration file... " << conffile << endl; 
  fconf.open( conffile, ios::in | ios::binary );
  fconf.read( (char*)&vlink[0][0], 8*nsite*sizeof(int) );
  fconf.read( (char*)&vllink[0][0], 8*nsite*sizeof(int) );
  fconf.read( (char*)&vplaq[0][0], 6*nsite*sizeof(int) );
  fconf.read( (char*)&vflux[0][0], 2*nsite*sizeof(int) );
  fconf.close( );
}

//_________________________________________________________________________
inline void meas_init_conf( )
{
  for ( int is=0; is<nsite; is++)
  {
    nflux[0][vflux[is][0]]++;
    nflux[1][vflux[is][1]]++;

    nlink[0] += vlink[is][3];
    nlink[1] += vlink[is][7];

    nplaq[abs(vplaq[is][0])]++;
    nplaq[abs(vplaq[is][1])]++;
    nplaq[abs(vplaq[is][2])]++;
    nplaq[abs(vplaq[is][3])]++;
    nplaq[abs(vplaq[is][4])]++;
    nplaq[abs(vplaq[is][5])]++;
  }
}

//_______________________________________________________________________
void init(int argc, char *argv[], int seed )
{
  read_params( argc, argv );
  rlxd_init( 1, seed );
  init_lattice( leng, vneib );

  memset ( nblock, 0, (3+3*LENGF)*sizeof(int) );
  
  if ( !rconf ) // cold start
  { 
    memset ( nplaq,  0, LENGF*sizeof(int) );
    memset ( &nflux[0][0],  0, 2*LENGF*sizeof(int) );
    memset ( nlink,  0, 2*sizeof(int) );

    nplaq[0] = 6*nsite;
    nflux[0][0] = nsite;
    nflux[1][0] = nsite;

    memset ( vlink,  0, nsite*8*sizeof(int) );
    memset ( vllink, 0, nsite*8*sizeof(int) );
    memset ( vplaq,  0, nsite*6*sizeof(int) );
    memset ( vflux,  0, nsite*2*sizeof(int) );
  }
  else // load configuration stored in conffile
  {
    read_conf( );
#ifdef CHECK
    check( nsite, vneib );
#endif
    meas_init_conf( );
  }

    // precompute the factorials
  for ( int i=0; i<LENGTH; i++ )
    fac[i] = log( fact( i ) );

}

#endif
