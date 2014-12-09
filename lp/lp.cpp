#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "problem_description.hpp"
#include "dictionary.hpp"
#include "bprinter/table_printer.h"

typedef bool(*test_func)(const std::string&);

static bool run_test1(const std::string& input_file)
{
    problem_description pd(input_file.c_str ());
    dictionary d(pd);
	unsigned t_enter, t_leave;
    const dictionary::result_t r = d.pivot(t_enter, t_leave);

    const std::string res = input_file + ".output";

    std::ifstream ifs(res);
    std::string first;
    ifs >> first;
    if (first[0] != 'U')
    {
        std::istringstream is(first);
        unsigned enter;
        is >> enter;
        unsigned leave;
        double z;
        ifs >> leave >> z;
        if (r == dictionary::UNBOUNDED)
        {
            std::cout << "FAILURE: " << input_file << std::endl;
            std::cout << "Expected: Bounded, actual: Unbounded";
            return false;
        }
        

        if (enter == t_enter && leave == t_leave && std::abs(z - d.z ()) < 0.001)
        {
            return true;
        }
        std::cout << "FAILURE: " << input_file << std::endl;
        bprinter::TablePrinter tb(&std::cout);
        tb.AddColumn("", 4);
        tb.AddColumn("Expected", 6);
        tb.AddColumn("Actual", 6);

        tb.PrintHeader();
        tb << "Entr" << enter << t_enter << bprinter::endl()
            << "Leav" << leave << t_leave << bprinter::endl()
            << "Objc" << z << d.z () << bprinter::endl();
        tb.PrintFooter();

        return false;
    }
    else
    {
        if (r != dictionary::UNBOUNDED)
        {
            std::cout << "Expected: Unbounded, actual: Bounded";
            return false;
        }
        return true;
    }
}

static bool run_test(const std::string& base_dir, int num, test_func f)
{
    std::ostringstream os;
    os << base_dir << "/dict" << num;

    return f(os.str());
}

static void do_unittest(const std::string& base_dir, test_func f)
{
    int good = 0;
    int bad = 0;
    for (int i = 1; i <= 10; i++)
    {
        if (run_test(base_dir, i, f))
        {
            good++;
        }
        else
        {
            bad++;
        }
    }

    std::cout << "Passed: " << good << "\tFAILED: " << bad << std::endl;
}

static void run_assignment1(const std::string& base_dir, int num)
{
    std::ostringstream os;
    os << base_dir << "/part" << num << ".dict";

    std::ostringstream of;
    of << base_dir << "/part" << num << ".output";

    std::ofstream ofs(of.str(), std::ios::trunc);

    problem_description pd(os.str ().c_str());
    dictionary d(pd);
	unsigned t_enter;
	unsigned t_leave;
    const dictionary::result_t r = d.pivot(t_enter, t_leave);
    std::cout << "ASSIGNMENT #" << num << std::endl;
    if (r == dictionary::UNBOUNDED)
    {
        ofs << "UNBOUNDED.";
        std::cout << "UNBOUNDED" << std::endl;
    }
    else
    {
		ofs << t_enter << "\n" << t_leave << "\n" << d.z();
		std::cout << t_enter << "\n" << t_leave << "\n" << d.z() << std::endl;
    }
    std::cout << "================================================" << std::endl;
}

// =================================================================================================
static bool run_test2(const std::string& input_file)
{
    problem_description pd(input_file.c_str());
	std::size_t num_steps;
	double dz;
    const bool r = dictionary::solve(pd, dz, &num_steps) == dictionary::eFINAL;

    const std::string res = input_file + ".output";

    std::ifstream ifs(res);
    std::string first;
    ifs >> first;
    if (first[0] != 'U')
    {
        std::istringstream is(first);
        std::size_t steps;
        
        double z;
        is >> z;
        ifs >> steps;
        if (r)
        {
            if (steps == num_steps && std::abs (dz - z) < 0.001)
            {
                return true;
            }
            std::cout << "FAILURE: " << input_file << std::endl;
            bprinter::TablePrinter tb(&std::cout);
            tb.AddColumn("", 4);
            tb.AddColumn("Expected", 6);
            tb.AddColumn("Actual", 6);

            tb.PrintHeader();
            tb << "Step" << steps << num_steps << bprinter::endl()
                << "Objc" << z << dz << bprinter::endl();
            tb.PrintFooter();
            return false;
        }
        else
        {
            std::cout << "FAILURE: " << input_file << std::endl;
            std::cout << "Expected bounded" << std::endl;
            return false;
        }
    }
    else
    {
        if (!r)
        {
            return true;
        }
        std::cout << "FAILURE: " << input_file << std::endl;
        std::cout << "Expected unbounded" << std::endl;
        return false;
    }
}

// -------------------------------------------------------------------------------
static void run_assignment2(const std::string& base_dir, int num)
{
    std::ostringstream os;
    os << base_dir << "/part" << num << ".dict";

    std::ostringstream of;
    of << base_dir << "/part" << num << ".output";

    std::ofstream ofs(of.str(), std::ios::trunc);

    problem_description pd(os.str().c_str());
	double result;
	std::size_t steps;
    const bool r = dictionary::solve(pd, result, &steps) == dictionary::eFINAL;
    std::cout << "ASSIGNMENT #" << num << std::endl;
    if (!r)
    {
        ofs << "UNBOUNDED.";
        std::cout << "UNBOUNDED" << std::endl;
    }
    else
    {
        ofs << result << std::endl << steps;
        std::cout << result << std::endl << steps << std::endl;
    }
    std::cout << "================================================" << std::endl;
}
// -------------------------------------------------------------------------------
bool run_test_3(const std::string& base_dir, int test_num)
{
	std::ostringstream os;
	os << base_dir << "/test" << test_num << ".dict";
	problem_description pd(os.str ().c_str ());
	double res;
	dictionary::solution_t rc = dictionary::solve(pd, res);

	std::ostringstream tos;
	tos << base_dir << "/test" << test_num << ".output";
	std::ifstream ifs(tos.str());
	std::string x;
	ifs >> x;

	dictionary::solution_t out_rc = dictionary::eFINAL;
	double out_res;
	if (x[0] == 'U')
	{
		out_rc = dictionary::eUNBOUNDED;
	}
	else
	{
		if (x[0] == 'I')
		{
			out_rc = dictionary::eINFEASIBLE;
		}
		else
		{
			std::istringstream is(x);
			is >> out_res;
		}
	} 
	if (out_rc != rc)
	{
		return false;
	}
	if (rc == dictionary::eFINAL)
	{
		double delta = std::abs(out_res - res);
		return delta < 0.2;
	}
	return true;
}
// -------------------------------------------------------------------------------
void run_assignment_3(const std::string& base_dir, int test_num)
{
	std::ostringstream os;
	os << base_dir << "/part" << test_num << ".dict";
	problem_description pd(os.str().c_str());
	double res;
	dictionary::solution_t rc = dictionary::solve(pd, res);
	std::ostringstream tos;
	tos << base_dir << "/part" << test_num << ".output";
	std::ofstream ofs(tos.str());

	if (rc == dictionary::eFINAL)
	{
		ofs << res;
	}
	else
	{
		if (rc == dictionary::eUNBOUNDED)
		{
			ofs << "UNBOUNDED";
		}
		else
		{
			ofs << "INFEASIBLE";
		}
	}

}
// -------------------------------------------------------------------------------
bool run_test_4 (const std::string& base_dir, int test_num, bool f)
{
	std::ostringstream os;
	if (!f)
	{
		os << base_dir << "/test" << test_num;
	}
	else
	{
		os << base_dir << "/ilpTest" << test_num;
	}
	problem_description pd(os.str().c_str());
	double res;
	dictionary::solution_t rc = dictionary::ilp_solve(pd, res);

	std::ostringstream tos;
	if (!f)
	{
		tos << base_dir << "/test" << test_num << ".output";
	}
	else
	{
		tos << base_dir << "/ilpTest" << test_num << ".output";
	}
	std::ifstream ifs(tos.str());
	std::string x;
	ifs >> x;

	dictionary::solution_t out_rc = dictionary::eFINAL;
	double out_res;
	if (x[0] == 'U' || x[0] == 'u')
	{
		out_rc = dictionary::eUNBOUNDED;
	}
	else
	{
		if (x[0] == 'I' || x[0] == 'u')
		{
			out_rc = dictionary::eINFEASIBLE;
		}
		else
		{
			std::istringstream is(x);
			is >> out_res;
		}
	}
	if (out_rc != rc)
	{
		return false;
	}
	if (rc == dictionary::eFINAL)
	{
		double delta = std::abs(out_res - res);
		return delta < 0.2;
	}
	return true;
}
// -------------------------------------------------------------------------------
void run_test_4(const char* basedir)
{
	for (int i = 0; i < 100; i++)
	{
		try
		{
			if (!run_test_4(basedir, i, false))
			{
				std::cout << "ERROR: " << i << std::endl;
			}
			else
			{
				std::cout << i << std::endl;
			}
		}
		catch (std::exception & e)
		{
			std::cout << i << ") " << e.what() << std::endl;
		}
	}
	for (int i = 0; i < 10; i++)
	{
		try
		{
			if (!run_test_4(basedir, i, true))
			{
				std::cout << "ERROR: ilp" << i << std::endl;
			}
			else
			{
				std::cout << "ilp" << i << std::endl;
			}
		}
		catch (std::exception & e)
		{
			std::cout << "ilp" << i << ") " << e.what() << std::endl;
		}
	}
	std::cout << "Done" << std::endl;
}
// -------------------------------------------------------------------------------
void run_assignment_4(const std::string& base_dir, int test_num)
{
	std::ostringstream os;
	os << base_dir << "/part" << test_num << ".dict";
	problem_description pd(os.str().c_str());
	double res;
	dictionary::solution_t rc = dictionary::ilp_solve(pd, res);
	std::ostringstream tos;
	tos << base_dir << "/part" << test_num << ".output";
	std::ofstream ofs(tos.str());

	if (rc == dictionary::eFINAL)
	{
		ofs << res;
	}
	else
	{
		if (rc == dictionary::eUNBOUNDED)
		{
			ofs << "UNBOUNDED";
		}
		else
		{
			ofs << "INFEASIBLE";
		}
	}

}
// ---------------------------------------------------------------------------
void run_assignment_4(const char* base)
{
	for (int i = 1; i <= 5; i++)
	{
		run_assignment_4(base, i);
	}
}
// -------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

	
	//run_test_4(argv[1]);
	run_assignment_4(argv[1]);
	
	return 0;
}
// -------------------------------------------------------------------------------
int main0(int argc, char* argv[])
{
	//run_test_3("d:\\proj\\lp\\part3TestCases\\unitTests\\50", 57);

	int failed[] = {100, 6, 10, 11, 13, 14, 18, 27, 31, 
					36, 41, 42, 43, 46, 48, 49, 60, 
					66, 68, 69, 73, 93, 96, 99
					,100 };
	for (int i = 1; i <= 10; i++)
	{
		std::cout << i;
		std::cout.flush();
		bool bad = false;
		for (int k = 0; k < 100; k++)
		{
			if (failed[k] == -1)
			{
				continue;
			}
			if (failed[k] == 100)
			{
				break;
			}
			if (failed[k] == i)
			{
				bad = true;
				break;
			}
		}
		if (!bad)
		{
			run_assignment_3("d:\\proj\\lp\\part3TestCases\\assignmentParts", i);
			 
		}
		else
		{
			std::cout << " FAIL" << std::endl;
		}
	}

    
    return 0;
}
