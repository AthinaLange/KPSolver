/// CSV File IO module.
/** @file
*/

#ifndef _CSV_H
#define _CSV_H

#include <string>
#include <vector>

/// Reads in data line-by-line from a csv file into a vector - assume file has no header 
/**
@param[in]	fn	Filename.
@param[out]	eta	Surface height \f$\eta\f$.
*/
void    read_wave_csv(

	const std::string&    fn,
        std::vector<double>&  eta);

/// Writes the free surface height \f$\eta\f$ to a csv file 
/**
@param[in]	fn	Filename.
@param[in]	nx	Dimension of meshgrid: x.
@param[in]	ny	Dimension of meshgrid: y.
@param[in]	eta	Surface height \f$\eta\f$.
*/
void    write_soliton_csv(

	const std::string&          fn,
        int                         nx,
        int                         ny,
        const std::vector<double>&  eta);

/// Writes the \f$\eta\f$ output to a csv file - can also write \f$\Delta x\f$ and \f$\Delta y\f$ to file 
/**
@param[in]	fn	Filename.
@param[in]	nx	Dimension of meshgrid: x.
@param[in]	ny	Dimension of meshgrid: y.
@param[in]	delx	Spacing in x.
@param[in]	dely	Spacing in y.
@param[in]	eta	Surface height \f$\eta\f$.
*/
void    write_wave_csv(

	const std::string&          fn,
        int                         nx,
        int                         ny,
        double                      delx,
        double                      dely,
        const std::vector<double>&  eta);

#endif
