/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/problemcomm.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// Class ProblemComunicator
//

#ifndef problemcommunicatormode_h
#define problemcommunicatormode_h

#ifdef __PARALLEL_MODE

/**
 * ProblemCommunicatorMode determines the valid mode.
 * The mode is used to set up communication pattern, which differ for
 * node and element cut algorithms.
 * Additional remote element mode hass been added to capture the case, when CommunicatorM is intended to
 * support remote element data exchange (for example when nonlocal material models are present).
 */
enum ProblemCommunicatorMode {
  ProblemCommMode__UNKNOWN_MODE, ProblemCommMode__NODE_CUT,
  ProblemCommMode__ELEMENT_CUT, ProblemCommMode__REMOTE_ELEMENT_MODE
};

#endif
#endif // problemcommunicatormode_h
