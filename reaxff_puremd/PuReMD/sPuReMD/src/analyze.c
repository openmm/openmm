/*----------------------------------------------------------------------
  SerialReax - Reax Force Field Simulator

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "analyze.h"

#include <ctype.h>

#include "box.h"
#include "list.h"
#include "tool_box.h"
#include "vector.h"


#define MAX_FRAGMENT_TYPES (100)


enum atoms
{
    O_ATOM = 0,
    N_ATOM = 1,
    SI_ATOM = 2,
};

enum molecule_type
{
    UNKNOWN = 0,
    WATER = 1,
};


typedef struct
{
    int atom_count;
    int atom_list[MAX_MOLECULE_SIZE];
    int mtypes[MAX_ATOM_TYPES];
} molecule;


/* copy bond list into old bond list */
static void Copy_Bond_List( reax_system *system, control_params *control,
        reax_list **lists )
{
    int i, j, top_old;
    reax_list *new_bonds = lists[BONDS];
    reax_list *old_bonds = lists[OLD_BONDS];

    for ( top_old = 0, i = 0; i < system->N; ++i )
    {
        Set_Start_Index( i, top_old, old_bonds );

        // fprintf( stdout, "%d: ", i );
        for ( j = Start_Index( i, new_bonds ); j < End_Index( i, new_bonds ); ++j )
            if ( new_bonds->bond_list[j].bo_data.BO >= control->bg_cut )
            {
                // fprintf( stderr, "%d ", new_bonds->bond_list[j].nbr );
                old_bonds->bond_list[ top_old ].nbr =
                    new_bonds->bond_list[j].nbr;
                old_bonds->bond_list[ top_old ].bo_data.BO =
                    new_bonds->bond_list[j].bo_data.BO;
                top_old++;
            }

        Set_End_Index( i, top_old, old_bonds);
        // fprintf( stderr, "--- s: %d, e: %d\n",
        // Start_Index( i, old_bonds ),  End_Index( i, old_bonds ) );
    }
}


// ASSUMPTION: Bond lists are sorted
static int Compare_Bond_Lists( int atom, control_params *control, reax_list **lists )
{
    int oldp, newp;
    reax_list *new_bonds = lists[BONDS];
    reax_list *old_bonds = lists[OLD_BONDS];

    /*fprintf( stdout, "\n%d\nold_bonds:", atom );
      for( oldp = Start_Index( atom, old_bonds );
           oldp < End_Index( atom, old_bonds ); ++oldp )
      if( old_bonds->bond_list[oldp].bo_data.BO >= control->bg_cut )
      fprintf( stdout, "%5d", old_bonds->bond_list[oldp].nbr );

      fprintf( stdout, "\nnew_bonds:" );
      for( newp = Start_Index( atom, new_bonds );
           newp < End_Index( atom, new_bonds ); ++newp )
      if( new_bonds->bond_list[newp].bo_data.BO >= control->bg_cut )
      fprintf( stdout, "%5d", new_bonds->bond_list[newp].nbr );*/


    for ( oldp = Start_Index( atom, old_bonds ),
            newp = Start_Index( atom, new_bonds );
            oldp < End_Index(atom, old_bonds) || newp < End_Index(atom, new_bonds);
            oldp = MIN( oldp + 1, End_Index( atom, old_bonds ) ),
            newp = MIN( newp + 1, End_Index( atom, new_bonds ) ) )
    {
        while ( oldp < End_Index( atom, old_bonds )
                && old_bonds->bond_list[oldp].bo_data.BO < control->bg_cut )
        {
            ++oldp;
        }

        while ( newp < End_Index( atom, new_bonds )
                && new_bonds->bond_list[newp].bo_data.BO < control->bg_cut )
        {
            ++newp;
        }

        /*fprintf( fout, "%d, oldp: %d - %d, newp: %d - %d",
          atom, oldp, old_bonds->bond_list[oldp].nbr,
          newp,  new_bonds->bond_list[newp].nbr );*/

        if ( oldp < End_Index( atom, old_bonds ) )
        {
            /* there are some other bonds in the old list */
            if ( newp < End_Index( atom, new_bonds ) )
            {
                if ( old_bonds->bond_list[oldp].nbr !=
                        new_bonds->bond_list[newp].nbr )
                {
                    //fprintf( fout, " --> case1, return 1\n" );
                    return 1;
                }
            }
            else
            {
                /* there is no other bond in the new list */
                //fprintf( fout, " --> case2, return 1\n" );
                return 1;
            }
        }
        else
        {
            /* there are no more bonds in old_bond list */
            if ( newp < End_Index( atom, new_bonds ) )
            {
                /* there is at least one other bond in the new list */
                //fprintf( fout, " --> case 3, return 1\n" );
                return 1;
            }
            else
            {
                /* there is no other bond in the new list, either */
                //fprintf( fout, " --> case 4, return 0\n" );
                return 0;
            }
        }
    }

    return 0;
}


static void Get_Molecule( int atom, molecule *m, int *mark, reax_system *system,
        control_params *control, reax_list *bonds, int print, FILE *fout )
{
    int i, start, end;

    start = Start_Index( atom, bonds );
    end = End_Index( atom, bonds );

    if ( print )
    {
        fprintf( fout, "%5d(%2s)",
                 atom + 1, system->reax_param.sbp[ system->atoms[atom].type ].name );
    }
    mark[atom] = 1;
    m->atom_list[ m->atom_count++ ] = atom;
    m->mtypes[ system->atoms[ atom ].type ]++;

    for ( i = start; i < end; ++i )
    {
        if ( bonds->bond_list[i].bo_data.BO >= control->bg_cut &&
                !mark[bonds->bond_list[i].nbr] )
        {
            Get_Molecule( bonds->bond_list[i].nbr, m, mark,
                          system, control, bonds, print, fout );
        }
    }
}


static void Print_Molecule( reax_system *system, molecule *m, int mode,
        char *s, size_t s_size, char *s_out, size_t s_out_size )
{
    int j, atom;

    s[0] = '\0';

    if ( mode == 1 )
    {
        /* print molecule summary */
        for ( j = 0; j < MAX_ATOM_TYPES; ++j )
        {
            if ( m->mtypes[j] )
            {
                snprintf( s_out, s_out_size, "%.*s%14s%3d", (int) MAX(s_out_size - 18, 0),
                        s, system->reax_param.sbp[j].name, m->mtypes[j] );
                s_out[s_out_size - 1] = '\0';
            }
        }
    }
    else if ( mode == 2 )
    {
        /* print molecule details */
        for ( j = 0; j < m->atom_count; ++j )
        {
            atom = m->atom_list[j];
            snprintf( s_out, s_out_size, "%.*s%.*s(%d)", (int) MAX((s_out_size - 5) / 2, 0), s,
                    (int) MAX((s_out_size - 5) / 2, 0),
                    system->reax_param.sbp[ system->atoms[atom].type ].name, atom );
            s_out[s_out_size - 1] = '\0';
        }
    }
}


static void Analyze_Molecules( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, FILE *fout )
{
    int atom, i, j, k, l, flag;
    int *mark = workspace->mark;
    int *old_mark = workspace->old_mark;
    int num_old, num_new;
    char s[MAX_MOLECULE_SIZE * 10];
    char s_out[MAX_MOLECULE_SIZE * 10];
    reax_list *new_bonds = lists[BONDS];
    reax_list *old_bonds = lists[OLD_BONDS];
    molecule old_molecules[20], new_molecules[20];

    fprintf( fout, "molecular analysis @ %d\n", data->step );
    for ( i = 0; i < system->N; i++ )
    {
        mark[i] = 0;
    }
    for ( i = 0; i < system->N; i++ )
    {
        old_mark[i] = 0;
    }

    /* compare new molecules to old molecules */
    for ( atom = 0; atom < system->N; ++atom )
        if ( !mark[atom] )
        {
            //fprintf( fout, "atom %d is being compared\n", atom );
            if ( Compare_Bond_Lists( atom, control, lists ) )
            {
                /* old & new lists are different, print the old and new molecules */
                num_old = num_new = 0;
                flag = 0;
                i = k = 0;
                memset( &old_molecules[0], 0, sizeof(molecule) * 20 );
                memset( &new_molecules[0], 0, sizeof(molecule) * 20 );

                fprintf( fout, "get_old_mol: " );
                Get_Molecule( atom, &old_molecules[num_old++], old_mark,
                              system, control, old_bonds, 1, fout );
                fprintf( fout, "\nget_new_mol: " );
                Get_Molecule( atom, &new_molecules[num_new++], mark,
                              system, control, new_bonds, 1, fout );
                fprintf( fout, "\n" );

                while ( !flag )
                {
                    flag = 1;

                    for ( ; i < num_old; ++i )
                        for ( j = 0; j < old_molecules[i].atom_count; ++j )
                            if ( !mark[old_molecules[i].atom_list[j]] )
                            {
                                fprintf( fout, "get_new_mol: " );
                                Get_Molecule( old_molecules[i].atom_list[j],
                                              &new_molecules[num_new++],
                                              mark, system, control, new_bonds, 1, fout );
                                fprintf( fout, "\n" );
                                flag = 0;
                            }

                    for ( ; k < num_new; ++k )
                        for ( l = 0; l < new_molecules[k].atom_count; ++l )
                            if ( !old_mark[new_molecules[k].atom_list[l]] )
                            {
                                fprintf( fout, "get_old_mol: " );
                                Get_Molecule( new_molecules[k].atom_list[l],
                                              &old_molecules[num_old++],
                                              old_mark, system, control, old_bonds, 1, fout );
                                fprintf( fout, "\n" );
                                flag = 0;
                            }
                }

                fprintf( fout, "old molecules: " );
                for ( i = 0; i < num_old; ++i )
                {
                    Print_Molecule( system, &old_molecules[i], 1, s,
                            sizeof(s), s_out, sizeof(s_out) );
                    fprintf( fout, "%s\t", s );
                }
                fprintf( fout, "\n" );

                fprintf( fout, "new molecules: " );
                for ( i = 0; i < num_new; ++i )
                {
                    Print_Molecule( system, &new_molecules[i], 1, s,
                            sizeof(s), s_out, sizeof(s_out) );
                    fprintf( fout, "%s\t", s );
                }
                fprintf( fout, "\n" );
            }
        }

    Copy_Bond_List( system, control, lists );

    //sfree( mark, "Analyze_Molecules::mark" );
    //sfree( old_mark, "Analyze_Molecules::old_mark" );

    fprintf( fout, "\n" );
    fflush( fout );
}


#if defined(DEBUG_FOCUS)
static void Report_Bond_Change( reax_system *system, control_params *control,
        static_storage *workspace,  reax_list *old_bonds,
        reax_list *new_bonds, int a1, int a2, int flag, FILE *fout )
{
    int i;
    int rev1, rev2;
    int mol1 = -1, mol2 = -1;
    // which molecule the atom belongs to, 0: Silica, 1: Water

    rev1 = workspace->orig_id[a1];
    rev2 = workspace->orig_id[a2];

    if ( !strncmp( system->atoms[a1].name, "  Si", 8 ) ||
            !strncmp( system->atoms[a1].name, "   O", 8 ) )
    {
        mol1 = 0;
    }
    else
    {
        mol1 = 1;
    }

    if ( !strncmp( system->atoms[a2].name, "  Si", 8 ) ||
            !strncmp( system->atoms[a2].name, "   O", 8 ) )
    {
        mol2 = 0;
    }
    else
    {
        mol2 = 1;
    }

    if ( mol1 == 0 && mol2 == 0 )   // silica-silica
    {
        if ( flag )
        {
            fprintf( fout, "silica bond formation:" );
        }
        else
        {
            fprintf( fout, "silica bond breakage :" );
        }

        fprintf( fout, "%5d(%s)-%5d(%s)\n",
                 rev1, system->atoms[a1].name, rev2, system->atoms[a2].name );
    }
    else if ( mol1 == 1 && mol2 == 1 )  // water-water
    {
        if ( flag )
        {
            fprintf( fout, "water bond formation:" );
        }
        else
        {
            fprintf( fout, "water bond breakage :" );
        }

        fprintf( fout, "%5d(%s)-%5d(%s)\n",
                 rev1, system->atoms[a1].name, rev2, system->atoms[a2].name );
    }
    else    // water-silica!
    {
        if ( flag )
        {
            fprintf( fout, "SILICA-WATER bond formation:" );
        }
        else
        {
            fprintf( fout, "SILICA-WATER bond breakage :" );
        }

        fprintf( fout, "%5d(%s)-%5d(%s)\n",
                 rev1, system->atoms[a1].name, rev2, system->atoms[a2].name );

        fprintf( fout, "%5d(%s) was connected to:", rev1, system->atoms[a1].name );
        for ( i = Start_Index(a1, old_bonds); i < End_Index(a1, old_bonds); ++i )
        {
            if ( old_bonds->bond_list[i].bo_data.BO >= control->bg_cut )
            {
                fprintf( fout, " %5d(%s)",
                         workspace->orig_id[ old_bonds->bond_list[i].nbr ],
                         system->atoms[ old_bonds->bond_list[i].nbr ].name );
            }
        }
        fprintf( fout, "\n" );

        fprintf( fout, "%5d(%s) was connected to:", rev2, system->atoms[a2].name );
        for ( i = Start_Index(a2, old_bonds); i < End_Index(a2, old_bonds); ++i )
        {
            if ( old_bonds->bond_list[i].bo_data.BO >= control->bg_cut )
            {
                fprintf( fout, " %5d(%s)",
                         workspace->orig_id[ old_bonds->bond_list[i].nbr ],
                         system->atoms[ old_bonds->bond_list[i].nbr ].name );
            }
        }
        fprintf( fout, "\n" );
    }
}
#endif


/* ASSUMPTION: Bond lists are sorted */
#if defined(DEBUG_FOCUS)
static void Compare_Bonding( int atom, reax_system *system, control_params *control,
        static_storage *workspace, reax_list *old_bonds,
        reax_list *new_bonds, FILE *fout )
{
    int oldp, newp;

    /* fprintf( fout, "\n%d\nold_bonds:", atom );
       for( oldp = Start_Index( atom, old_bonds );
            oldp < End_Index( atom, old_bonds ); ++oldp )
       if( old_bonds->bond_list[oldp].bo_data.BO >= control->bg_cut )
       fprintf( fout, "%5d", old_bonds->bond_list[oldp].nbr );

       fprintf( fout, "\nnew_bonds:" );
       for( newp = Start_Index( atom, new_bonds );
            newp < End_Index( atom, new_bonds ); ++newp )
       if( new_bonds->bond_list[newp].bo_data.BO >= control->bg_cut )
       fprintf( fout, "%6d", new_bonds->bond_list[newp].nbr );
       fprintf( fout, "\n" ); */

    for ( oldp = Start_Index( atom, old_bonds );
            oldp < End_Index( atom, old_bonds )
            && old_bonds->bond_list[oldp].nbr < atom; ++oldp )
        ;

    for ( newp = Start_Index( atom, new_bonds );
            newp < End_Index( atom, new_bonds )
            && new_bonds->bond_list[newp].nbr < atom; ++newp )
        ;

    while ( oldp < End_Index( atom, old_bonds ) ||
            newp < End_Index( atom, new_bonds ) )
    {
        while ( oldp < End_Index( atom, old_bonds ) &&
                old_bonds->bond_list[oldp].bo_data.BO < control->bg_cut )
        {
            ++oldp;
        }

        while ( newp < End_Index( atom, new_bonds ) &&
                new_bonds->bond_list[newp].bo_data.BO < control->bg_cut )
        {
            ++newp;
        }

        /*fprintf( fout, "%d, oldp: %d - %d: %f    newp: %d - %d: %f",
          atom, oldp, old_bonds->bond_list[oldp].nbr,
          old_bonds->bond_list[oldp].bo_data.BO,
          newp,  new_bonds->bond_list[newp].nbr,
          new_bonds->bond_list[newp].bo_data.BO ); */

        if ( oldp < End_Index( atom, old_bonds ) )
        {
            /* there are some more bonds in the old list */
            if ( newp < End_Index( atom, new_bonds ) )
            {
                if ( old_bonds->bond_list[oldp].nbr <
                        new_bonds->bond_list[newp].nbr )
                {
                    // fprintf( fout, "%5d-%5d bond broken\n",
                    // atom, old_bonds->bond_list[oldp].nbr );
                    Report_Bond_Change( system, control, workspace, old_bonds, new_bonds,
                                        atom, old_bonds->bond_list[oldp].nbr, 0,
                                        fout );
                    ++oldp;
                }
                else if ( old_bonds->bond_list[oldp].nbr >
                          new_bonds->bond_list[newp].nbr )
                {
                    // fprintf( fout, "%5d-%5d bond formed\n",
                    // atom, new_bonds->bond_list[newp].nbr );
                    Report_Bond_Change( system, control, workspace, old_bonds, new_bonds,
                                        atom, new_bonds->bond_list[newp].nbr, 1,
                                        fout );
                    ++newp;
                }
                else
                {
                    ++newp;
                    ++oldp;
                }
            }
            else
                /* there is no other bond in the new list */
                while ( oldp < End_Index( atom, old_bonds ) )
                {
                    if ( old_bonds->bond_list[oldp].bo_data.BO >= control->bg_cut )
                    {
                        // fprintf( fout, "%5d-%5d bond broken\n",
                        // atom, old_bonds->bond_list[oldp].nbr );
                        Report_Bond_Change( system, control, workspace,
                                            old_bonds, new_bonds, atom,
                                            old_bonds->bond_list[oldp].nbr, 0,
                                            fout );
                    }
                    ++oldp;
                }
        }
        else
        {
            /* there are no more bonds in old_bond list */
            if ( newp < End_Index( atom, new_bonds ) )
            {
                /* there is at least one other bond in the new list */
                while ( newp < End_Index( atom, new_bonds ) )
                {
                    if ( new_bonds->bond_list[newp].bo_data.BO >= control->bg_cut )
                    {
                        // fprintf( fout, "%5d-%5d bond formed\n",
                        // atom, new_bonds->bond_list[newp].nbr );
                        Report_Bond_Change( system, control, workspace,
                                old_bonds, new_bonds, atom,
                                new_bonds->bond_list[newp].nbr, 1, fout );
                    }
                    ++newp;
                }
            }
            else
            {
                /* there is no other bond in the new list, either --
                   no need to do anything */
            }
        }
    }
}
#endif


static void Visit_Bonds( int atom, int *mark, int *type, reax_system *system,
                  control_params *control, reax_list *bonds, int ignore )
{
    int i, t, start, end, nbr;
    real bo;

    mark[atom] = 1;
    t = system->atoms[atom].type;
    if ( ignore && control->ignore[t] )
    {
        return;
    }
    type[t]++;

    start = Start_Index( atom, bonds );
    end = End_Index( atom, bonds );
    for ( i = start; i < end; ++i )
    {
        nbr = bonds->bond_list[i].nbr;
        bo = bonds->bond_list[i].bo_data.BO;
        if ( bo >= control->bg_cut && !mark[nbr] )
            Visit_Bonds( nbr, mark, type, system, control, bonds, ignore );
    }
}


static void Analyze_Fragments( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, FILE *fout, int ignore )
{
    int atom, i, flag;
    int *mark = workspace->mark;
    int num_fragments, num_fragment_types;
    char fragment[MAX_ATOM_TYPES];
    char fragment_out[MAX_ATOM_TYPES];
    char fragments[MAX_FRAGMENT_TYPES][MAX_ATOM_TYPES];
    int fragment_count[MAX_FRAGMENT_TYPES];
    molecule m;
    reax_list *new_bonds = lists[BONDS];
//    reax_list *old_bonds = lists[OLD_BONDS];

    /* fragment analysis */
    fprintf( fout, "step%d fragments\n", data->step );
    num_fragments = 0;
    num_fragment_types = 0;
    for ( i = 0; i < system->N; i++ )
    {
        mark[i] = 0;
    }

    for ( atom = 0; atom < system->N; ++atom )
    {
        if ( !mark[atom] )
        {
            /* discover a new fragment */
            memset( m.mtypes, 0, MAX_ATOM_TYPES * sizeof(int) );
            Visit_Bonds( atom, mark, m.mtypes, system, control, new_bonds, ignore );
            ++num_fragments;
            Print_Molecule( system, &m, 1, fragment, sizeof(fragment),
                    fragment_out, sizeof(fragment_out) );

            /* check if a similar fragment already exists */
            flag = 0;
            for ( i = 0; i < num_fragment_types; ++i )
            {
                if ( !strncmp( fragments[i], fragment_out, MAX_ATOM_TYPES ) )
                {
                    ++fragment_count[i];
                    flag = 1;
                    break;
                }
            }

            if ( flag == 0 )
            {
                /* it is a new one, add to the fragments list */
                strncpy( fragments[num_fragment_types], fragment_out,
                        sizeof(fragments[num_fragment_types]) );
                fragments[num_fragment_types][sizeof(fragments[num_fragment_types]) - 1] = '\0';
                fragment_count[num_fragment_types] = 1;
                ++num_fragment_types;
            }
        }
    }

    /* output the results of fragment analysis */
    for ( i = 0; i < num_fragment_types; ++i )
    {
        if ( strnlen(fragments[i], MAX_ATOM_TYPES) )
        {
            fprintf( fout, "%d of %s\n", fragment_count[i], fragments[i] );
        }
    }
    fprintf( fout, "\n" );
    fflush( fout );

    /* compare new bonds to old bonds */
    //for( atom = 0; atom < system->N; ++atom ) {
    // fprintf( fout, "atom: %d\n", atom ); fflush( fout );
    // Compare_Bonding( atom, system, control, workspace,
    // old_bonds, new_bonds, fout );
    //}
    //Copy_Bond_List( system, control, lists );
}


#if defined(DEBUG_FOCUS)
static void Analyze_Silica( reax_system *system, control_params *control,
                     simulation_data *data, static_storage *workspace,
                     reax_list **lists, FILE *fout )
{
    int atom, i, j, k, pi, pk, pk_j, newp, coord;
    int O_SI_O_count, SI_O_SI_count;
    int si_coord[10], ox_coord[10];
    real O_SI_O, SI_O_SI;
    reax_list *new_bonds = lists[BONDS];
    reax_list *thb_intrs = lists[THREE_BODIES];

    Analyze_Fragments( system, control, data, workspace, lists, fout, 0 );

    /* analyze atom coordinations */
    for ( i = 0; i < 10; ++i )
    {
        si_coord[i] = ox_coord[i] = 0;
    }

    for ( atom = 0; atom < system->N; ++atom )
    {
        coord = 0;

        for ( newp = Start_Index( atom, new_bonds );
                newp < End_Index( atom, new_bonds ); ++newp )
        {
            if ( new_bonds->bond_list[newp].bo_data.BO >= control->bg_cut )
            {
                ++coord;
            }
        }

        if ( system->atoms[ atom ].type == SI_ATOM )
        {
            /*if( coord == 4 )
            full_coord_SI++;
            else less_coord_SI++;*/
            ++si_coord[coord];
        }
        else if ( system->atoms[ atom ].type == O_ATOM )
        {
            /*if( coord == 2 )
            full_coord_O++;
            else less_coord_O++;*/
            ++ox_coord[coord];
        }
    }

    /* fprintf( fout, "\nFour Coordinated SI: %.2f%\n",
       (double)full_coord_SI / (full_coord_SI + less_coord_SI) * 100. );
       fprintf( fout, "Four Coordinated O : %.2f%\n",
       (double)full_coord_O  / (full_coord_O  + less_coord_O ) * 100. ); */

    fprintf( fout, "Silicon coordinations:\n" );
    for ( i = 1; i < 10; ++i )
        if ( si_coord[i] )
            fprintf( fout, "\t%d-coord: %d\n", i, si_coord[i] );

    fprintf( fout, "\nOxygen coordinations:\n" );
    for ( i = 1; i < 10; ++i )
        if ( ox_coord[i] )
            fprintf( fout, "\t%d-coord: %d\n", i, ox_coord[i] );


    /* analyze bond angles */
    O_SI_O = 0;
    O_SI_O_count = 0;

    SI_O_SI = 0;
    SI_O_SI_count = 0;

    for ( j = 0; j < system->N; ++j )
    {
        if ( system->atoms[j].type == O_ATOM || system->atoms[j].type == SI_ATOM )
        {
            for ( pi = Start_Index(j, new_bonds); pi < End_Index(j, new_bonds); ++pi )
            {
                if ( new_bonds->bond_list[pi].bo_data.BO >= control->bg_cut )
                {
                    i = new_bonds->bond_list[pi].nbr;

                    if (system->atoms[i].type == O_ATOM || system->atoms[i].type == SI_ATOM)
                    {
                        for ( pk = Start_Index( pi, thb_intrs );
                                pk < End_Index( pi, thb_intrs ); ++pk )
                        {
                            k = thb_intrs->three_body_list[pk].thb;
                            pk_j = thb_intrs->three_body_list[pk].pthb;
                            // get k's pointer on j's bond list

                            if ( new_bonds->bond_list[pk_j].bo_data.BO >=
                                    control->bg_cut )   // physical j&k bond
                            {
                                /*fprintf( fout, "%5d(%d) %5d(%d) %5d(%d)   %8.3f\n",
                                  i, system->atoms[i].type, j, system->atoms[j].type,
                                  k, system->atoms[k].type,
                                  thb_intrs->three_body_list[pk].theta );*/

                                if ( system->atoms[i].type == O_ATOM &&
                                        system->atoms[j].type == SI_ATOM &&
                                        system->atoms[k].type == O_ATOM )
                                {
                                    O_SI_O_count++;
                                    O_SI_O += thb_intrs->three_body_list[pk].theta;
                                }
                                else if ( system->atoms[i].type == SI_ATOM &&
                                          system->atoms[j].type == O_ATOM &&
                                          system->atoms[k].type == SI_ATOM )
                                {
                                    SI_O_SI_count++;
                                    SI_O_SI += thb_intrs->three_body_list[pk].theta;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fprintf( fout, "\nAverage O-Si-O angle: %8.2f\n",
             RAD2DEG(O_SI_O / O_SI_O_count) );
    fprintf( fout, "Average Si-O-Si angle:  %8.2f\n\n\n",
             RAD2DEG(SI_O_SI / SI_O_SI_count) );

    fflush( fout );
}
#endif


static int Get_Type_of_Molecule( molecule *m )
{
    if ( m->atom_count == 3 && m->mtypes[1] == 2 && m->mtypes[2] == 1 )
    {
        return WATER;
    }

    return UNKNOWN;
}


static void Calculate_Dipole_Moment( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list *bonds, FILE *fout )
{
    int i, atom, count;
    molecule m;
    real mu_sum;
    rvec tmpvec, mu;
    int *mark = workspace->mark;

    mu_sum = 0;
    count = 0;

    for ( atom = 0; atom < system->N; ++atom )
    {
        /* start discovering water molecules from the central O atom */
        if ( !mark[atom] && system->atoms[atom].type == 2 )
        {
            rvec_MakeZero( mu );

            memset( &m, 0, sizeof(molecule) );
            Get_Molecule( atom, &m, mark, system, control, bonds, 0, fout );

            if ( Get_Type_of_Molecule( &m ) == WATER )
            {
                ++count;

                for ( i = 1; i < 2; ++i )
                {
                    Compute_Atom_Distance_Triclinic( control,
                            &system->box,
                            system->atoms[ m.atom_list[0] ].x,
                            system->atoms[ m.atom_list[i] ].x,
                            system->atoms[ m.atom_list[0] ].rel_map,
                            system->atoms[ m.atom_list[i] ].rel_map,
                            system->atoms[ m.atom_list[0] ].rel_map, //TODO: what to use for rel_box?
                            tmpvec );
                    rvec_ScaledAdd( mu, -system->atoms[m.atom_list[0]].q / 2.0, tmpvec );
                }

                mu_sum += rvec_Norm( mu );
            }
        }
    }

    fprintf( fout, "%7d  %10d      %10.5f\n",
             data->step, count, mu_sum / count * ECxA_to_DEBYE );
    fflush( fout );
}


static void Copy_Positions( reax_system *system, static_storage *workspace )
{
    int i;

    for ( i = 0; i < system->N; ++i )
    {
        rvec_Copy( workspace->x_old[i], system->atoms[i].x );
    }
}


static void Calculate_Drift( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, FILE *fout )
{
    int i, type;
    int count[MAX_ATOM_TYPES];
    real drift;
    rvec driftvec;
    real sum_sqr_drift[MAX_ATOM_TYPES], sum_drift[MAX_ATOM_TYPES];

    for ( i = 0; i < MAX_ATOM_TYPES; ++i )
    {
        count[i] = 0;
        sum_sqr_drift[i] = sum_drift[i] = 0.0;
    }

    for ( i = 0; i < system->N; ++i )
    {
        //if( control->restrict_type == -1 ||
        // system->atoms[i].type == control->restrict_type )
        if ( workspace->x_old[i][0] > -system->box.box_norms[0] &&
                workspace->x_old[i][1] > -system->box.box_norms[1] &&
                workspace->x_old[i][2] > -system->box.box_norms[2] )
        {
            type = system->atoms[i].type;
            ++count[type];

            Compute_Atom_Distance_Triclinic( control, &system->box,
                    workspace->x_old[i], system->atoms[i].x,
                    system->atoms[i].rel_map, system->atoms[i].rel_map,
                    system->atoms[i].rel_map, driftvec ); //TODO: what to use for rel_box?

            if ( FABS( driftvec[0] ) >= system->box.box_norms[0] / 2.0 - 2.0 ||
                    FABS( driftvec[1] ) >= system->box.box_norms[0] / 2.0 - 2.0 ||
                    FABS( driftvec[2] ) >= system->box.box_norms[0] / 2.0 - 2.0 )
            {
                /* the atom has moved almost half the box size.
                   exclude it from further drift computations as it might have an
                   improper contribution due to periodic boudnaries. */
                workspace->x_old[i][0] = -999999999999.0;
                workspace->x_old[i][1] = -999999999999.0;
                workspace->x_old[i][2] = -999999999999.0;
                continue;
            }

            drift = rvec_Norm_Sqr( driftvec );
            sum_sqr_drift[type] += drift;
            sum_drift[type]     += SQRT( drift );
        }
    }

    fprintf( fout, "%7d  oxy  %6d  %10.6f\n",
             data->step, count[2], sum_sqr_drift[2] / count[2] );
    fprintf( fout, "%7d  hyd  %6d  %10.6f\n",
             data->step, count[1], sum_sqr_drift[1] / count[1] );

    fflush( fout );
}


#if defined(DEBUG_FOCUS)
static void Calculate_Density_3DMesh( reax_system *system, simulation_data *data,
        FILE *fout )
{
    int i, j, k;
    int occupied_cells;
    int ***cell_counter;
    ivec my_cell;
    rvec mesh_cell_lens = { 1, 1, 1 };
    ivec mesh_dims;
    real occupied_vol, density;

    /* determine the mesh dimensions based on the current box size */
    for ( i = 0; i < 3; ++i )
    {
        mesh_dims[i] = system->box.box_norms[i] / mesh_cell_lens[i] + 0.99;
    }

    fprintf( stderr, "mesh_dims: %3d  %3d  %3d\n",
             mesh_dims[0], mesh_dims[1], mesh_dims[2] );

    /* allocate counter for each mesh cell */
    cell_counter = (int ***) scalloc( mesh_dims[0], sizeof(int),
           "Calculate_Density_3DMesh::cell_counter" );

    for ( i = 0; i < mesh_dims[0]; ++i )
    {
        cell_counter[i] = (int **) scalloc( mesh_dims[1], sizeof(int),
               "Calculate_Density_3DMesh::cell_counter[i]" );

        for ( j = 0; j < mesh_dims[1]; ++j )
        {
            cell_counter[i][j] = (int *) scalloc( mesh_dims[2], sizeof(int),
                   "Calculate_Density_3DMesh::cell_counter[i][j]" );
        }
    }


    /* go over the atoms and assign each into their corresponding mesh cell */
    for ( i = 0; i < system->N; ++i )
    {
        my_cell[0] = system->atoms[i].x[0] / mesh_cell_lens[0];
        my_cell[1] = system->atoms[i].x[1] / mesh_cell_lens[1];
        my_cell[2] = system->atoms[i].x[2] / mesh_cell_lens[2];

        cell_counter[ my_cell[0] ][ my_cell[1] ][ my_cell[2] ]++;
    }


    /* calculate volume occupied */
    occupied_cells = 0;
    for ( i = 0; i < mesh_dims[0]; ++i )
    {
        for ( j = 0; j < mesh_dims[1]; ++j )
        {
            for ( k = 0; k < mesh_dims[2]; ++k )
            {
                if ( cell_counter[i][j][k] )
                {
                    ++occupied_cells;
                }
            }
        }
    }

    occupied_vol = occupied_cells * mesh_cell_lens[0]
        * mesh_cell_lens[1] * mesh_cell_lens[2];

    fprintf( stderr, "occupied cells: %8d\n", occupied_cells );
    fprintf( stderr, "occupied vol  : %8.2f\n", occupied_vol );
    fprintf( stderr, "system volume : %8.2f\n", system->box.volume );
    fprintf( stderr, "system mass   : %8.2f\n", data->M );

    density = data->M * AMU_to_GRAM / (occupied_vol * POW( ANG_to_CM, 3.0 ));
    fprintf( stderr, "density       : %g\n", density );
    fprintf( stderr, "AMU_to_GRAM   : %g\n", AMU_to_GRAM );
    fprintf( stderr, "ANG_to_CM     : %g\n", ANG_to_CM );
}
#endif


#if defined(DEBUG_FOCUS)
static void Calculate_Density_Slice( reax_system *system, simulation_data *data, FILE *fout )
{
    real slice_thickness = 0.5;
    int *slice_occ;
    int i, num_slices, my_slice, max_occ = 0;

    /* allocate counter */
    num_slices = system->box.box_norms[2] / slice_thickness + 1.;
    slice_occ = (int*) scalloc( num_slices, sizeof(int),
           "Calculate_Density_Slice::slice_occ" );

    /* distribute atoms to slices */
    for ( i = 0; i < system->N; ++i )
    {
        my_slice = system->atoms[i].x[2] / slice_thickness;
        slice_occ[ my_slice ]++;
    }

    /* determine maximum occupancy */
    for ( i = 0; i < num_slices; ++i )
    {
        fprintf( stderr, "occ[%d]: %d\n", i, slice_occ[i] );
        if ( slice_occ[i] > max_occ )
        {
            max_occ = slice_occ[i];
        }
    }

    /* find luzzatti-interface slice */
    for ( i = 0; i < num_slices; ++i )
    {
        if ( (real)slice_occ[i] / max_occ > 0.5 )
        {
            fprintf( stderr, "%d - %d is the luzzatti interface\n", i - 1, i );
            break;
        }
    }
}
#endif


void Analysis( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int steps;

    steps = data->step - data->prev_steps;

    /****** Molecular Analysis ******/
    if ( control->molec_anal
            && steps % control->freq_molec_anal == 0 )
    {
        if ( steps == 1 )
        {
            if ( lists[OLD_BONDS]->allocated == FALSE )
            {
                Make_List( lists[BONDS]->n, lists[BONDS]->n_max, lists[BONDS]->total_intrs,
                        TYP_BOND, lists[OLD_BONDS] );
            }

            if ( control->molec_anal == REACTIONS )
            {
                Copy_Bond_List( system, control, lists );
            }
            if ( control->diffusion_coef )
            {
                Copy_Positions( system, workspace );
            }
        }

        if ( control->molec_anal == FRAGMENTS )
        {
            /* discover molecules */
            Analyze_Fragments( system, control, data, workspace, lists,
                    out_control->mol, 0 );
            /* discover fragments without the ignored atoms */
            if ( control->num_ignored )
            {
                Analyze_Fragments( system, control, data, workspace, lists,
                        out_control->ign, 1 );
            }
        }
        else if ( control->molec_anal == REACTIONS )
        {
            /* discover molecular changes - reactions */
            Analyze_Molecules( system, control, data, workspace, lists,
                    out_control->mol );
        }
    }

    /****** Electric Dipole Moment ******/
    if ( control->dipole_anal && steps % control->freq_dipole_anal == 0 )
    {
        Calculate_Dipole_Moment( system, control, data, workspace,
                lists[BONDS], out_control->dpl );
    }

    /****** Drift ******/
    if ( control->diffusion_coef && steps % control->freq_diffusion_coef == 0 )
    {
        Calculate_Drift( system, control, data, workspace, out_control->drft );
    }
}
