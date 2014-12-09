#ifndef __DICTIONARY_HPP__
#define __DICTIONARY_HPP__

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "problem_description.hpp"

class dictionary
{
public:
    enum result_t
    {
        DONE,
        CONT,
        UNBOUNDED
    };

	enum solution_t 
	{
		eFINAL,
		eUNBOUNDED,
		eINFEASIBLE
	};
public:
    dictionary(const problem_description& pd);

    void print(std::ostream& os) const;

    bool feasible() const;
    double b(unsigned r) const;
    double z() const;

	static solution_t solve(const problem_description& pd, double& result, std::size_t* steps = 0);
	static solution_t solve(dictionary& d, double& result, const problem_description& pd, std::size_t* steps = 0);

	static solution_t ilp_solve(const problem_description& pd, double& result);

    result_t pivot(unsigned& var_enter, unsigned& var_leave);
	result_t pivot();

    std::size_t steps() const;
	void dualize();
private:
	typedef std::vector <double> row_t;
private:
    void find_entering_variable(std::vector <unsigned>& candidates) const;
    // returns saturation value (if any)
    double find_leaving_variable_indecies(unsigned column, std::vector <unsigned>& indecies) const;

    struct leave_info
    {
        unsigned      basic_variable;
        std::size_t   non_basic_variable;
        double        saturation_value;
    };
    // returns index of the first unbounded result, or info.size ()
    std::size_t find_leaving_variables_set(std::vector <unsigned>& entering, std::vector <leave_info>& info) const;

    // returns index in the info vector
    std::size_t choose_entering_and_leaving_variable(const std::vector <leave_info>& info) const;

    void exchange(const leave_info& inf);

	std::string get_basic_var(unsigned i) const;
	std::string get_non_basic_var(unsigned i) const;

	void change_basis(const std::vector <unsigned>& B, const std::vector <unsigned>& I, const std::vector <double>& orig_coeffs);

	void restore_original_objective(const row_t& original_coeffs, const problem_description& pd);
	void get_fractional_rows(std::vector <unsigned>& rows) const;
	void add_cutting_planes(const std::vector <unsigned>& rows);
private:
    unsigned m;
    unsigned n;
    std::vector <row_t> matrix;
    std::vector <unsigned> basic;
    std::vector <unsigned> non_basic;
    row_t coeffs;
    
	std::size_t steps_num;

	std::vector <unsigned> complimentary_mapping;
	bool is_primal;
};

// ======================================================================================

inline 
std::ostream& operator << (std::ostream& os, const dictionary& d)
{
    d.print(os);
    return os;
}

inline
std::ostream& operator << (std::ostream& os, dictionary::result_t r)
{
    switch (r)
    {
    case dictionary::DONE:
        os << "DONE";
        break;
    case dictionary::CONT:
        os << "CONT";
        break;
	case dictionary::UNBOUNDED:
        os << "UNBOUNDED";
        break;
    }
    return os;
}

#endif