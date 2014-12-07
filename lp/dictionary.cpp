#include <iostream>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <algorithm>

#include "dictionary.hpp"
#include "bprinter/table_printer.h"
#include "matrix2.hxx"



//#define DEBUG

//#define DEBUG_INI

static const double EPSILON = 0.00001;
// -----------------------------------------------------------------
static bool is_zero(double x)
{
    
    return std::abs(x) < EPSILON;
}
// -----------------------------------------------------------------
static bool eq(double x, double y)
{
    return is_zero(x - y);
}
// -----------------------------------------------------------------
static bool gt(double x, double y)
{
    if (eq(x, y))
    {
        return false;
    }
    return x > y;
}
// -----------------------------------------------------------------
static bool gte(double x, double y)
{
    return  (eq(x, y) || (x > y));
}
// -----------------------------------------------------------------
static bool lt(double x, double y)
{
    return !gte(x, y);
}
// -----------------------------------------------------------------
static bool lte(double x, double y)
{
    return eq(x, y) || lt(x, y);
}
// -----------------------------------------------------------------
dictionary::dictionary(const problem_description& pd)
    : m(pd.m),
    n(pd.n),
    basic(pd.basic),
    non_basic(pd.non_basic),
	is_primal(true),
	steps_num(0)
{
    matrix.reserve(m);
    unsigned k = 0;
    for (unsigned r = 0; r < m; r++)
    {
        row_t row;
        row.reserve(n + 1);
        row.push_back(pd.b[r]);
        for (unsigned c = 0; c < n; c++)
        {
            row.push_back(pd.a[k++]);
        }
        matrix.push_back(row);
    }
    coeffs.reserve(n + 1);
    coeffs.push_back(pd.z0);
    for (unsigned c = 0; c < n; c++)
    {
        coeffs.push_back(pd.c[c]);
    }
}
// -------------------------------------------------------------------
double dictionary::b(unsigned r) const
{
    return matrix[r][0];
}
// -------------------------------------------------------------------
double dictionary::z() const
{
    return coeffs[0];
}
// -------------------------------------------------------------------
bool dictionary::feasible() const
{
    bool ok = true;
    for (unsigned r = 0; r < m; r++)
    {
        if (!gte(b(r), 0.0))
        {
            ok = false;
            break;
        }
    }
    return ok;
}
// -------------------------------------------------------------------
void dictionary::find_entering_variable(std::vector <unsigned>& candidates) const
{
    for (unsigned c = 1; c < n + 1; c++)
    {
        if (gt(coeffs[c], 0))
        {
            candidates.push_back(c - 1);
        }
    }
}
// -------------------------------------------------------------------
double dictionary::find_leaving_variable_indecies(unsigned column, std::vector <unsigned>& indecies) const
{
    unsigned index = m;
    double minimum = -1.0;
    bool first_time = true;
    bool is_unbounded = true;
    for (unsigned r = 0; r < m; r++)
    {
        double br = b(r);
        double a = matrix[r][column];

        if (gte(br, 0.0) && lt(a, 0.0))
        {
            double k = -br / a;
            if (first_time)
            {
                minimum = k;
                index = r;
                first_time = false;
            }
            else
            {
                if (lt(k, minimum))
                {
                    minimum = k;
                    index = r;
                }
            }
        }
    }
    if (index == m)
    {
        return 0.0;
    }
    for (unsigned r = 0; r < m; r++)
    {
        double br = b(r);
        double a = matrix[r][column];

        if (gte(br, 0.0) && lt(a, 0.0))
        {
            double k = -br / a;
            if (eq(k, minimum))
            {
                indecies.push_back(r);
            }
        }
    }
    return minimum;
}
// -------------------------------------------------------------------
std::size_t dictionary::find_leaving_variables_set(std::vector <unsigned>& entering, std::vector <leave_info>& info) const
{
    std::size_t unbounded_idx = entering.size();
    
    std::vector <unsigned> basic_indecies;

    for (std::size_t c = 0; c < entering.size(); c++)
    {
        const unsigned column = entering[c] + 1;
        assert(column < n + 1);
        const double saturation = find_leaving_variable_indecies(column, basic_indecies);
        if (basic_indecies.empty())
        {
            unbounded_idx = c;
            break;
        }
        for (std::size_t j = 0; j < basic_indecies.size(); j++)
        {
            leave_info li;
            li.basic_variable = basic_indecies[j];
            li.non_basic_variable = entering[c];
            li.saturation_value = saturation;
            info.push_back(li);
        }
        basic_indecies.clear();
    }



    return unbounded_idx;
}
// -------------------------------------------------------------------
std::size_t dictionary::choose_entering_and_leaving_variable(const std::vector <leave_info>& info) const
{
    // blands rule
    
    assert(!info.empty());

    unsigned min_non_basic_name = non_basic [info [0].non_basic_variable];
    std::size_t min_non_basic_index = 0;
    
    for (std::size_t i = 1; i < info.size(); i++)
    {
        const unsigned name = non_basic[info[i].non_basic_variable];
        if (name < min_non_basic_name)
        {
            min_non_basic_name = name;
            min_non_basic_index = i;
        }
    }
    
    unsigned min_basic_name = basic[info[min_non_basic_index].basic_variable];
    std::size_t min_basic_index = min_non_basic_index;

    for (std::size_t i = 1; i < info.size(); i++)
    {
        const unsigned non_basic_name = non_basic[info[i].non_basic_variable];
        if (non_basic_name == min_non_basic_name)
        {
            const unsigned basic_name = basic[info[i].basic_variable];
            if (basic_name < min_basic_name)
            {
                min_basic_name = basic_name;
                min_basic_index = i;
            }
        }
    }

    return min_basic_index;
}
// -------------------------------------------------------------------
void dictionary::exchange(const leave_info& inf)
{
    row_t new_row(n + 1);
    row_t& old_row = matrix[inf.basic_variable];
    unsigned column = inf.non_basic_variable + 1;
    for (unsigned c = 0; c < n + 1; c++)
    {
        new_row[c] = -old_row[c] / old_row[column];
    }
	new_row[column] = 1.0 / old_row[column];

    matrix[inf.basic_variable] = new_row;
    for (unsigned r = 0; r < m; r++)
    {
        if (r == inf.basic_variable)
        {
            continue;
        }
        double k = matrix[r][column];
        if (is_zero(k))
        {
            continue;
        }
        for (unsigned c = 0; c < n + 1; c++)
        {
            if (c != column)
            {
                matrix[r][c] = matrix[r][c] + k*new_row[c];
            }
            else
            {
                matrix[r][c] = k*new_row[c];
            }
        }
    }
    double k = coeffs[column];
    for (unsigned c = 0; c < n + 1; c++)
    {
        if (c != column)
        {
            coeffs[c] = coeffs[c] + k*new_row[c];
        }
        else
        {
            coeffs[c] = k*new_row[c];
        }
    }
    std::swap(basic[inf.basic_variable], non_basic[inf.non_basic_variable]);
}
// -------------------------------------------------------------------
static inline 
std:: size_t index_of(const std::vector <unsigned>& array, unsigned value)
{
	std::vector <unsigned>::const_iterator itr = std::lower_bound(array.begin(), array.end(), value);
	assert(itr != array.end());
	assert(*itr == value);
	return itr - array.begin();
}
// -------------------------------------------------------------------
typedef techsoft::matrix<double> matrix_t;
void check_row(const matrix_t& m, unsigned r)
{
	assert(r < m.rowno());
}
void check_col(const matrix_t& m, unsigned c)
{
	assert(c < m.colno());
}
// -------------------------------------------------------------------
void dictionary::change_basis(const std::vector <unsigned>& B, const std::vector <unsigned>& I, const std::vector <double>& old_coeffs)
{
	// create problem matrix
	
	matrix_t A(m, m + n);
	for (unsigned r = 0; r < m; r++)
	{
		check_row(A, r);

		techsoft::Mat_iter<double> row = A.row(r);
		for (unsigned c = 0; c < n; c++)
		{
			check_col(A, c);
			row[c] = matrix[r][c + 1];//pd.a[pd.m*r + c];
		}
	}
	for (unsigned r = 0; r < m; r++)
	{
		check_row(A, r);
		techsoft::Mat_iter<double> row = A.row(r);

		check_col(A, n+r);
		row[n + r] = -1.0;
	}

	matrix_t Bv(m, 1);
	for (unsigned r = 0; r < m; r++)
	{
		check_row(Bv, r);
		check_col(Bv, 0);
		Bv.row(r)[0] = matrix[r][0];
	}

	typedef std::map <unsigned, unsigned> columns_mapping_t;
	columns_mapping_t columns_mapping;

	for (std::size_t i = 0; i < n; i++)
	{
		columns_mapping[non_basic[i]] = i;
	}

	for (std::size_t i = 0; i < m; i++)
	{
		columns_mapping[basic[i]] = n + i;
	}

	matrix_t C(m + n, 1);
	
	
	for (unsigned j = 0; j < n; j++)
	{
		check_row(C, j);
		check_col(C, 0);
		C.row(j)[0] = old_coeffs[j];
		
	}
	
	
	matrix_t AB(m, m);
	matrix_t AI(m, n);
	matrix_t CB(1, m);
	for (std::size_t i = 0; i < m; i++)
	{
		unsigned col = columns_mapping[B[i]];
		check_row(C, col);
		check_col(C, 0);

		double v = C.row(col)[0];
		check_row(CB, 0);
		check_col(CB, i);
		CB.row(0)[i] = v;
		for (std::size_t r = 0; r < m; r++)
		{
			check_row(A, r);
			check_col(A, col);

			check_row(AB, r);
			check_col(AB, i);
			AB.column(i)[r] = A.column(col)[r];
		}
	}
	matrix_t CI(1, I.size());
	for (std::size_t i = 0; i < n; i++)
	{
		unsigned col = columns_mapping[I[i]];
		check_row(CI, 0);
		check_col(CI, i);

		check_row(C, col);
		check_col(C, 0);

		CI.row(0)[i] = C.row(col)[0];
		for (std::size_t r = 0; r < m; r++)
		{
			check_row(A, r);
			check_col(A, col);

			check_row(AI, r);
			check_col(AI, i);
			AI.column(i)[r] = A.column(col)[r];
		}
	}

	

	bool rc = AB.inv();
	assert(rc == true);

	matrix_t T = -AB*AI;
	matrix_t cm = CI + CB*T;
	matrix_t U = AB*Bv;
	matrix_t Z = CB*U;

	for (unsigned r = 0; r < m; r++)
	{
		for (unsigned c = 1; c <= n; c++)
		{
			check_row(T, r);
			check_col(T, c-1);

			matrix[r][c] = T.row(r)[c - 1];
		}
		check_row(U, r);
		check_col(U, 0);
		matrix[r][0] = U.row(r)[0];
	}
	
	double dz = Z.row(0)[0];
	coeffs[0] = dz;
	for (std::size_t i = 0; i < I.size(); i++)
	{
		check_row(cm, 0);
		check_col(cm, i);
		coeffs[i + 1] = cm.row(0)[i];
	}
	basic = B;
	non_basic = I;
	m = basic.size();
	n = non_basic.size();
}


void dictionary::restore_original_objective(const row_t& original_coeffs, const problem_description& pd)
{
	std::vector <double> new_coeffs;
	std::vector <unsigned> B = basic;
	std::vector <unsigned> I = non_basic;

	
	for (int i = 0; i < pd.n; i++)
	{
		new_coeffs.push_back(-1);
	}
	change_basis(pd.basic, pd.non_basic, new_coeffs);
	
	
	change_basis(B, I, pd.c);

	
}
// -------------------------------------------------------------------
dictionary::solution_t dictionary::solve(const problem_description& pd, double& result, std::size_t* steps)
{
	dictionary d(pd);

	if (!d.feasible ())
	{
#if defined(DEBUG_INI)
		std::cout << "PRIMAL" << std::endl;
		std::cout << d << std::endl;
#endif

		row_t original_coeffs = d.coeffs;
		std::vector <unsigned> original_non_basic = d.non_basic;
		for (unsigned i = 1; i <= d.n; i++)
		{
			d.coeffs[i] = -1.0;
		}
		d.dualize ();
#if defined(DEBUG_INI)		
		std::cout << "DUAL" << std::endl;
		std::cout << d << std::endl;
#endif
		assert(d.feasible());
		
		while (true)
		{
			const result_t rc = d.pivot();
			if (rc != CONT)
			{
				if (rc == UNBOUNDED)
				{
					return eINFEASIBLE;
				}
				else
				{
					break;
				}
			}
		}
#if defined(DEBUG_INI)		
		std::cout << "Final dual" << std::endl;
		std::cout << d << std::endl;
#endif
		
		d.dualize ();
#if defined(DEBUG_INI)		
		std::cout << "Primal with biased objective" << std::endl;
		std::cout << d << std::endl;
#endif
		assert(d.feasible());
		d.restore_original_objective(original_coeffs, pd);
#if defined(DEBUG_INI)		
		std::cout << "Primal with resotred objective" << std::endl;
		std::cout << d << std::endl;
#endif
	}

    while (true)
    {
#if defined(DEBUG)
        std::cout << d << std::endl;
        std::cout << "Feasible: " << d.feasible() << std::endl;
#endif
        const result_t rc = d.pivot();
        if (rc != CONT)
        {
			if (rc == DONE)
			{
				if (steps)
				{
					*steps = d.steps();
				}
				result = d.z();
				return eFINAL;
			}
			else
			{
				return eUNBOUNDED;
			}
        }
    }
	assert(false);
	return eFINAL;
}
// -------------------------------------------------------------------
dictionary::result_t dictionary::pivot()
{
	unsigned x;
	unsigned y;
	return pivot(x, y);
}
// -------------------------------------------------------------------
dictionary::result_t dictionary::pivot(unsigned& var_enter, unsigned& var_leave)
{
    std::vector <unsigned> entering;
    find_entering_variable(entering);
#if defined(DEBUG)

    std::cout << *this << std::endl;

    if (!entering.empty())
    {
        std::cout << "Entering candidates: ";
        for (std::size_t i = 0; i < entering.size(); i++)
        {
            std::cout << get_non_basic_var (entering[i]) << " ";
        }
        std::cout << std::endl;
    }
#endif

    if (entering.empty())
    {
#ifdef DEBUG
        std::cout << "No entering candidates found. Finished" << std::endl;
        std::cout << "Objective function: " << z() << std::endl;
        
#endif
        return DONE;
    }


    std::vector <leave_info> leaving;
    const std::size_t unbounded_idx = find_leaving_variables_set(entering, leaving);

    if (unbounded_idx < entering.size())
    {
#if defined(DEBUG)
        std::cout << "Dictionary is unbounded at " << get_non_basic_var (entering[unbounded_idx]) << std::endl;
#endif
        return UNBOUNDED;
    }

    const std::size_t pair_idx = choose_entering_and_leaving_variable(leaving);
    assert(pair_idx < leaving.size());

	const std::size_t non_bsc = leaving[pair_idx].non_basic_variable;
	const std::size_t bsc = leaving[pair_idx].basic_variable;
#if defined(DEBUG)
    std::cout << "Entering/Leaving pairs:" << std::endl;
    for (std::size_t c = 0; c < leaving.size(); c++)
    {
        if (c == pair_idx)
        {
            std::cout << "*\t";
        }
        else
        {
            std::cout << " \t";
        }
		const std::size_t cnon_bsc = leaving[c].non_basic_variable;
		const std::size_t cbsc = leaving[c].basic_variable;    
        std::cout << get_non_basic_var (cnon_bsc)  << " enters, "
            << get_basic_var (cbsc) << " leaves. Saturation: " << leaving[c].saturation_value << std::endl;
    }
#endif
	var_enter = non_basic[non_bsc];
	var_leave = basic[bsc];
	exchange(leaving[pair_idx]);
	steps_num++;
    return CONT;
}
// -------------------------------------------------------------------
std::size_t dictionary::steps() const
{
	return steps_num;
}
// -------------------------------------------------------------------
void dictionary::dualize()
{
	std::vector <row_t> tr(n);
	for (std::size_t r = 0; r < n; r++)
	{
		tr[r].resize(m + 1);
	}

	for (std::size_t r = 0; r < m; r++)
	{
		for (std::size_t c = 0; c < n; c++)
		{
			tr[c][r+1] = -matrix[r][c+1];
		}
	}

	for (std::size_t i = 0; i < n; i++)
	{
		tr[i][0] = -coeffs[i + 1];
	}
	double z0 = -coeffs[0];
	coeffs.resize(m + 1);
	coeffs[0] = z0;
	for (std::size_t i = 0; i < m; i++)
	{
		coeffs[i+1] = -matrix[i][0];
	}
	complimentary_mapping.resize(m + n);
	for (std::size_t i = 0; i < n; i++)
	{
		complimentary_mapping[i] = non_basic[i];
	}
	for (std::size_t i = 0; i < m; i++)
	{
		complimentary_mapping[i+n] = basic[i];
	}
	std::swap(basic, non_basic);
	std::swap(m, n);
	matrix = tr;
	is_primal ^= true;
}
// ---------------------------------------------------------------------------
std::string dictionary::get_basic_var(unsigned i) const
{
	std::ostringstream os;

	if (is_primal)
	{
		os << "X" << basic[i];
	}
	else
	{
		os << "Y";
		std::size_t p = 0;
		for (std::size_t j = 0; j < m + n; j++)
		{
			if (complimentary_mapping[j] == basic[i])
			{
				p = j + 1;
				break;
			}
		}
		assert(p != 0);
		os << p;
	}
	return os.str();
}
// ---------------------------------------------------------------------------
std::string dictionary::get_non_basic_var(unsigned i) const
{
	std::ostringstream os;

	if (is_primal)
	{
		os << "X" << non_basic[i];
	}
	else
	{
		os << "Y";
		std::size_t p = 0;
		for (std::size_t j = 0; j < m + n; j++)
		{
			if (complimentary_mapping[j] == non_basic[i])
			{
				p = j + 1;
				break;
			}
		}
		assert(p != 0);
		os << p;
	}
	return os.str();
}
// ---------------------------------------------------------------------------
void dictionary::print(std::ostream& os) const
{
    bprinter::TablePrinter tb(&os);
    tb.AddColumn("", 8);
    for (unsigned i = 0; i < n+1; i++)
    {

        tb.AddColumn("", 8);
    }
    tb.PrintFooter();
    for (unsigned r = 0; r < m; r++)
    {
        std::ostringstream os;
		os << get_basic_var(r);
        tb << os.str();
        os.str("");
        os.clear();
        const double b = matrix[r][0];
        if (!is_zero(b))
        {
            tb << b;
        }
        else
        {
            tb << "";
        }
        for (unsigned c = 1; c < n + 1; c++)
        {
            double a = matrix[r][c];
            if (!is_zero(a))
            {
                if (eq(a, 1.0))
                {
					os << get_non_basic_var(c - 1);
                }
                else
                {
                    if (eq(a, -1.0))
                    {
						os << "-" << get_non_basic_var(c - 1);
                    }
                    else
                    {
						os << a << get_non_basic_var(c - 1);
                    }
                }
                tb << os.str();
                os.str("");
                os.clear();
            }
            else
            {
                tb << "";
            }
        }
       
    }
    tb << bprinter::endl();
    tb.PrintFooter();
    tb << "z" << coeffs[0];
    for (unsigned c = 1; c < n + 1; c++)
    {
        double a = coeffs[c];
        if (!is_zero(a))
        {
            std::ostringstream os;
            if (eq(a, 1.0))
            {
				os << get_non_basic_var(c - 1);
            }
            else
            {
                if (eq(a, -1.0))
                {
					os << "-" << get_non_basic_var(c - 1);
                }
                else
                {
					os << a << get_non_basic_var(c - 1);
                }
            }
            tb << os.str();

        }
    }
    tb << bprinter::endl ();
    tb.PrintFooter();
}
