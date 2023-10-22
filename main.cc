#include <iomanip>
#include <filesystem>

#include "HyWall.h"
#include "PTL.h"

std::vector<std::string> Split(const std::string& str, char c)
{
	std::vector<std::string> output;
	std::string elem = "";
	for (int i = 0; i < str.length(); i++)
	{
		char cur = str[i];
		if (cur==c)
		{
			output.push_back(elem);
			elem = "";
		}
		else
		{
			elem += cur;
		}
	}
	output.push_back(elem);
	return output;
}

std::string cleandata(std::vector<std::string>& args)
{
	std::string output;
	bool foundnum = false;
	std::string content = args[0];
	std::string acceptables = "1234567890";
	for (int i = 0; i < content.length(); i++)
	{
		foundnum = foundnum || (acceptables.find(content[i]) != std::string::npos);
		if (foundnum) output += content[i];
	}
	return output;
}

std::string specialfread(std::vector<std::string>& args)
{
	int col_num = 0;
	std::istringstream iss(args[2]);
	iss >> col_num;
	std::vector<std::string> args2;
	args2.push_back(args[0]);
	args2.push_back(args[1]);
	std::string line = PTL::BuiltIns::PTLFunc_fread(args2);
	auto v = Split(line, ',');
	return v[col_num];
}

std::string suth(std::vector<std::string>& args)
{
	double tin = PTL::BuiltIns::AssertConvertDouble(args[0]);
	const double tref = 110.4;
	const double muref = 1.45151376745308e-06;
	double result = muref*pow(tin,1.5)/(tin+tref);
	return std::to_string(result);
}

struct input_data_t
{
	double dist, x, y, z, rho, mu, p, u, v, w, T, mu_t, mom_rhs, dp_dx;
};

struct output_data_t
{
	double tau, qw, vort, fail;
};

struct alt_settings_t
{
	bool compute_derivatives, compute_commutation;
	PTL::PropertySection* section;
	double eps_T, eps_U, eps_R;
	std::string commutation_filename;

	void read(PTL::PropertySection& section_in)
	{
		section = &section_in;
		section_in["compute_derivatives"].MapTo(&compute_derivatives) = new PTL::PTLBoolean(false, "Compute sensitivity hessian");
		section_in["compute_commutation"].MapTo(&compute_commutation) = new PTL::PTLBoolean(false, "Compute commutation error (requires DNS data)");
		
		section_in["commutation_filename"].MapTo(&commutation_filename) = new PTL::PTLString("qf.dat", "file name for dns data for commutation error");
		
		section_in["eps_T"].MapTo(&eps_T) = new PTL::PTLDouble(1e-4, "Finite difference perturbation for temperature");
		section_in["eps_U"].MapTo(&eps_U) = new PTL::PTLDouble(1e-4, "Finite difference perturbation for velocity");
		section_in["eps_R"].MapTo(&eps_R) = new PTL::PTLDouble(1e-4, "Finite difference perturbation for density");
	}
	void parse()
	{
		section->StrictParse();
	}
};

bool streq(const std::string& s1, const std::string& s2) {return s1.compare(s2)==0;}
void ReadInput(
	const std::string& filename,
	decltype(HyWall::settings)& settings,
	input_data_t& in_data,
	bool& ran_init,
	std::string& out_file,
	int argc,
	char** argv,
	alt_settings_t& alt_settings)
{
	PTL::AddUserDefinedFunction("cleandata", cleandata);
	PTL::AddUserDefinedFunction("specialfread", specialfread);
	PTL::AddUserDefinedFunction("suth", suth);
	PTL::PropertyTree input;
	PTL::Interactive i(argc, argv, &input);
	i.Run();
	ran_init = false;
	bool failed = false;
	bool is_init = false;
	std::string failed_message;
	try
	{
		input.Read(filename);
	}
	catch (PTL::PTLException e)
	{
		failed = true;
		is_init = filename == "--init";
		if (!is_init)
		{
			failed_message = e.what();
			std::cout << "FAILURE MESSAGE:\n" << failed_message << std::endl;
		}
	}
	
	settings.enableWallModel = true;
	settings.readRestart = false;
	input["WallModel"]["rayDim"].MapTo(&settings.rayDim)                   = new PTL::PTLInteger(30, "number of ray points");
	settings.asyncSolve = false;
	input["WallModel"]["verboseLevel"].MapTo(&settings.verboseLevel)          = new PTL::PTLInteger(1, "debug output level");
	input["WallModel"]["maxIterations"].MapTo(&settings.maxIterations)        = new PTL::PTLInteger(100, "Max. iterations");
	input["WallModel"]["wallSpacing"].MapTo(&settings.wallSpacing)            = new PTL::PTLDouble(1e-6, "Max. iterations");
	input["WallModel"]["wallTemperature"].MapTo(&settings.wallTemperature)    = new PTL::PTLDouble(100, "Wall Temperature");
	input["WallModel"]["adiabaticWall"].MapTo(&settings.adiabaticWall)        = new PTL::PTLBoolean(true, "Adiabatic wall");
	input["WallModel"]["variablePrandtlT"].MapTo(&settings.variablePrandtlT)  = new PTL::PTLBoolean(false, "Variable turbulent prandtl number");
	input["WallModel"]["fluidCp"].MapTo(&settings.fluidCp)                    = new PTL::PTLDouble(1005.0, "Specific heat");
	input["WallModel"]["turbPradntl"].MapTo(&settings.turbPradntl)            = new PTL::PTLDouble(0.72, "Turbulent Prandtl");
	input["WallModel"]["fluidPrandtl"].MapTo(&settings.fluidPrandtl)          = new PTL::PTLDouble(0.9, "Laminar Prandtl");
	input["WallModel"]["vanDriestAPlus"].MapTo(&settings.vanDriestAPlus)      = new PTL::PTLDouble(17.0, "van Driest Constant");
	input["WallModel"]["gasConstant"].MapTo(&settings.gasConstant)            = new PTL::PTLDouble(287.0, "Gas constant");
	input["WallModel"]["enableTransitionSensor"].MapTo(&settings.enableTransitionSensor) = new PTL::PTLBoolean(false, "Enable Transition Sensor");
	
	std::string mom_eq_str, trb_eq_str, eng_eq_str;
	input["WallModel"]["momentumEquationType"].MapTo(&mom_eq_str)          = new PTL::PTLString("ODE", "Momentum equation type");
	input["WallModel"]["turbulenceEquationType"].MapTo(&trb_eq_str)        = new PTL::PTLString("vanDriest", "Turbulence equation type");
	input["WallModel"]["energyEquationType"].MapTo(&eng_eq_str)            = new PTL::PTLString("ODE", "Energy equation type");
	
	input["WallModel"]["momentumUnderRelaxationODE"].MapTo(&settings.momentumUnderRelaxationODE)     = new PTL::PTLDouble(0.8, "Relaxation factor for momentum ODE");
	input["WallModel"]["turbulenceUnderRelaxationODE"].MapTo(&settings.turbulenceUnderRelaxationODE) = new PTL::PTLDouble(0.6, "Relaxation factor for turbulence ODE");
	input["WallModel"]["energyUnderRelaxationODE"].MapTo(&settings.energyUnderRelaxationODE)         = new PTL::PTLDouble(0.7, "Relaxation factor for energy ODE");
	input["WallModel"]["includeMomentumRhs"].MapTo(&settings.includeMomentumRhs)                     = new PTL::PTLBoolean(false, "Include the parameterized convection term");
	input["WallModel"]["isCompressible"].MapTo(&settings.isCompressible)                             = new PTL::PTLBoolean(false, "Use variable density");
	input["WallModel"]["suthViscRef"].MapTo(&settings.suthViscRef)                                   = new PTL::PTLDouble(1.45151376745308e-06, "Reference viscosity for viscosity power law");
	input["WallModel"]["suthTRef"].MapTo(&settings.suthTRef)                                         = new PTL::PTLDouble(110.4, "Reference temperature for viscosity power law");
	input["WallModel"]["errorTolerance"].MapTo(&settings.errorTolerance)                             = new PTL::PTLDouble(1e-6, "Error tolerance");
	
	std::string visc_law_str;
	input["WallModel"]["viscousLaw"].MapTo(&visc_law_str)                                            = new PTL::PTLString("constant", "Viscous law");
	
	input["InputData"]["dist"].MapTo(&in_data.dist)       = new PTL::PTLDouble(0.01,     "Sampling distance from wall");
	input["InputData"]["x"].MapTo(&in_data.x)             = new PTL::PTLDouble(0.0,      "X-coordinate of sampling location");
	input["InputData"]["y"].MapTo(&in_data.y)             = new PTL::PTLDouble(0.0,      "Y-coordinate of sampling location");
	input["InputData"]["z"].MapTo(&in_data.z)             = new PTL::PTLDouble(0.0,      "Z-coordinate of sampling location");
	input["InputData"]["rho"].MapTo(&in_data.rho)         = new PTL::PTLDouble(1.17683,  "Density at sampling location");
	input["InputData"]["mu"].MapTo(&in_data.mu)           = new PTL::PTLDouble(1e-5,     "Viscosity at sampling location");
	input["InputData"]["P"].MapTo(&in_data.p)             = new PTL::PTLDouble(101325.0, "Pressure at sampling location");
	input["InputData"]["U"].MapTo(&in_data.u)             = new PTL::PTLDouble(69.54,    "U-velocity at sampling location");
	input["InputData"]["V"].MapTo(&in_data.v)             = new PTL::PTLDouble(0.05,     "V-velocity at sampling location");
	input["InputData"]["W"].MapTo(&in_data.w)             = new PTL::PTLDouble(0.0,      "W-velocity at sampling location");
	input["InputData"]["T"].MapTo(&in_data.T)             = new PTL::PTLDouble(300.0,    "Temperature at sampling location");
	input["InputData"]["mu_t"].MapTo(&in_data.mu_t)       = new PTL::PTLDouble(0.0,      "Turbulent viscosity at sampling location");
	input["InputData"]["mom_rhs"].MapTo(&in_data.mom_rhs) = new PTL::PTLDouble(0.0,      "Momentum RHS at sampling location");
	input["InputData"]["dp_dx"].MapTo(&in_data.dp_dx)     = new PTL::PTLDouble(0.0,      "Pressure gradient at sampling location");
	
	std::string yscale_str;
	input["WallModel"]["yScale"].MapTo(&yscale_str) = new PTL::PTLString("trettelLarsson", "y-coordinate scaling");
	
	input["Output"]["filename"].MapTo(&out_file)    = new PTL::PTLString("solution.csv", "Output filename");
	
	alt_settings.read(input["Misc"]);
	
	if (!failed)
	{
		i.Run();
		input.StrictParse();
		alt_settings.parse();
	}
	     if (streq(mom_eq_str, "allmaras")) {settings.momentumEquationType = HyCore::momentum::allmaras;}
	else if (streq(mom_eq_str, "ODE"))      {settings.momentumEquationType = HyCore::momentum::ODE;}
	else if (streq(mom_eq_str, "fromFile")) {settings.momentumEquationType = HyCore::momentum::fromFile;}
	else if (!failed) {std::cout << "Invalid value for momentum equation" << std::endl; abort();}

	     if (streq(trb_eq_str, "linear"))    {settings.turbulenceEquationType = HyCore::turbulence::linear;}
	else if (streq(trb_eq_str, "ODE"))       {settings.turbulenceEquationType = HyCore::turbulence::ODE;}
	else if (streq(trb_eq_str, "vanDriest")) {settings.turbulenceEquationType = HyCore::turbulence::vanDriest;}
	else if (streq(trb_eq_str, "fromFile"))  {settings.turbulenceEquationType  = HyCore::turbulence::fromFile;}
	else if (streq(trb_eq_str, "pnlm"))      {settings.turbulenceEquationType  = HyCore::turbulence::pnlm;}
	else if (!failed) {std::cout << "Invalid value \"" << trb_eq_str << "\"for turbulence equation" << std::endl; abort();}

	     if (streq(eng_eq_str, "croccoBusemann")) {settings.energyEquationType = HyCore::energy::croccoBusemann;}
	else if (streq(eng_eq_str, "ODE"))            {settings.energyEquationType = HyCore::energy::ODE; }
	else if (streq(eng_eq_str, "linear"))         {settings.energyEquationType = HyCore::energy::linear;}
	else if (streq(eng_eq_str, "fromFile"))       {settings.energyEquationType = HyCore::energy::fromFile;}
	else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}

		if  (streq(yscale_str, "trettelLarsson")) {settings.yscaleType = HyCore::yscale::trettelLarsson;}
	else if (streq(yscale_str, "yPlus"))          {settings.yscaleType = HyCore::yscale::yPlus;}
	else if (streq(yscale_str, "mixed"))          {settings.yscaleType = HyCore::yscale::mixed;}
	else if (!failed) {std::cout << "Invalid value for y-coordinate transform" << std::endl; abort();}

	     if (streq(visc_law_str, "constant"))   {settings.viscousLaw = HyCore::visclaw::constant;}
	else if (streq(visc_law_str, "sutherland")) {settings.viscousLaw = HyCore::visclaw::sutherland;}
	else if (streq(visc_law_str, "PowerLaw"))   {settings.viscousLaw = HyCore::visclaw::PowerLaw;}
	else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}
	
	else if (is_init)
	{
		ran_init = true;
		std::cout << "Initializing input file..." << std::endl;
		input.CreateDefaultValuesFile("input.ptl");
		std::cout << "See \"input.ptl\"" << std::endl;
		exit(0);
	}
	else
	{
		std::cout << "Cannot find input file \"" << filename << "\". Try running with first argument \"--init\" instead." << std::endl;
		abort();
	}
}

struct hy_probes_t
{
	int u_prb, T_prb, mu_t_prb, y_prb, rho_prb, mu_prb, lam_t_prb;
};

struct hy_derivs_t
{
	
};

template <typename data_t> auto compute_deriv(const int du, const int dt, const int dr, const data_t& data, const double eu, const double et, const double er)
{
	std::vector<double> d0{0.0,  0.0,  1.0,  0.0,  0.0};
	std::vector<double> d1{0.0, -0.5,  0.0,  0.5,  0.0};
	std::vector<double> d2{0.0,  1.0, -2.0,  1.0,  0.0};
	std::vector<double> d3{0.5, -1.0,  0.0,  1.0, -0.5};
	std::vector<double> d4{1.0, -4.0,  6.0, -4.0,  1.0};
	
	std::vector<std::vector<double>> coeffs{d0, d1, d2, d3, d4};
	const double norm_u   = 1.0/std::pow(eu, du);
	const double norm_t   = 1.0/std::pow(et, dt);
	const double norm_r   = 1.0/std::pow(er, dr);
	const double norm_tot = norm_u*norm_t*norm_r;
	
	const auto coeffu = coeffs[du];
	const auto coeffr = coeffs[dr];
	const auto coefft = coeffs[dt];
	
	const int npt = HyWall::settings.rayDim;
	std::vector<double> deriv_E(npt, 0.0);
	std::vector<double> deriv_M(npt, 0.0);
	
	int nderiv = d4.size()-1;
	for (int ir = 0; ir <= nderiv; ++ir)
	{
		for (int it = 0; it <= nderiv; ++it)
		{
			for (int iu = 0; iu <= nderiv; ++iu)
			{
				std::vector<double> local_E(npt, 0.0);
				std::vector<double> local_M(npt, 0.0);
				const std::vector<double>& lam_t = std::get<6>(data[ir][it][iu]);
				const std::vector<double>& mu_t  = std::get<5>(data[ir][it][iu]);
				const std::vector<double>& mu    = std::get<4>(data[ir][it][iu]);
				const std::vector<double>& u     = std::get<1>(data[ir][it][iu]);
				const std::vector<double>& y     = std::get<0>(data[ir][it][iu]);
				const std::vector<double>& T     = std::get<2>(data[ir][it][iu]);
				for (int i = 0; i < npt; ++i)
				{
					int il = i-1;
					int iu = i+1;
					if (i==0)     il = i;
					if (i==npt-1) iu = i;
					
					double yl   = y[il];
					double yu   = y[iu];
					double uu   = u[iu];
					double ul   = u[il];
					double Tu   = T[iu];
					double Tl   = T[il];
					double mu_total_h   = 0.5*(mu[il]+mu[iu]) + 0.5*(mu_t[il]+mu_t[iu]);
					double lam_total_h  = HyWall::settings.fluidCp*(0.5*(mu[il]+mu[iu])/HyWall::settings.fluidPrandtl) + 0.5*(lam_t[il]+lam_t[iu]);
					double dudy = (uu-ul)/(yu-yl);
					double dTdy = (Tu-Tl)/(yu-yl);
					
					double lM = lam_total_h*dTdy;
					double lE = mu_total_h*dudy;
					local_E[i] = lE;
					local_M[i] = lM;
				}
				
				const double coeffloc = coeffr[ir]*coefft[it]*coeffu[iu]*norm_tot;
				
				for (int i = 0; i < npt; ++i)
				{
					deriv_M[i] += coeffloc*local_M[i];
					deriv_E[i] += coeffloc*local_E[i];
				}
			}
		}
	}
	return std::make_tuple(deriv_M, deriv_E);
}

void compute_wm_derivatives(hy_probes_t& probes, alt_settings_t& alt_settings, std::vector<double>& q_f)
{
	int tmp =  HyWall::settings.verboseLevel;
	HyWall::settings.verboseLevel = 0;
	std::cout << "computing derivatives..." << std::endl;
	using vcd = std::vector<double>;
	//mu, mu_t, lam_t, u, T
	auto perturb = [&](double dU, double dT, double dR, double T_f) -> void
	{
		q_f[0] += dR*T_f*HyWall::settings.gasConstant; //pressure
		q_f[1] += dU;
		q_f[4] += dT;
	};
	using hy_solution_t = std::tuple<vcd, vcd, vcd, vcd, vcd, vcd, vcd>;
	auto solve = [&](double dU, double dT, double dR, double T_f)-> hy_solution_t
	{
		perturb(dU, dT, dR, T_f);
		//continue here
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		std::vector<double> v_y      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_mu     (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_mu_t   (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_lam_t  (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_u      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_T      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_rho    (HyWall::settings.rayDim, 0.0);
		
		HyWall::ProbeSolution(probes.y_prb,    0, &v_y    [0]);
		HyWall::ProbeSolution(probes.u_prb,    0, &v_u    [0]);
		HyWall::ProbeSolution(probes.T_prb,    0, &v_T    [0]);
		HyWall::ProbeSolution(probes.rho_prb,  0, &v_rho  [0]);
		HyWall::ProbeSolution(probes.mu_prb,   0, &v_mu   [0]);
		HyWall::ProbeSolution(probes.mu_t_prb, 0, &v_mu_t [0]);
		HyWall::ProbeSolution(probes.lam_t_prb,0, &v_lam_t[0]);
		perturb(-dU, -dT, -dR, T_f);
		return std::make_tuple(v_y, v_u, v_T, v_rho, v_mu, v_mu_t, v_lam_t);
	};
	const int i_y     = 0;
	const int i_u     = 0;
	const int i_rho   = 0;
	const int i_T     = 0;
	const int i_mu    = 0;
	const int i_mu_t  = 0;
	const int i_lam_t = 0;
	
	const double P_f   = q_f[0];
	const double U_f   = q_f[1];
	const double T_f   = q_f[4];
	const double rho_f = P_f/(HyWall::settings.gasConstant*T_f);
	const double e_u = alt_settings.eps_U*U_f;
	const double e_t = alt_settings.eps_T*T_f;
	const double e_r = alt_settings.eps_R*rho_f;
	
	const int ngrid = 5;
	std::vector<hy_solution_t> sols_0(ngrid);
	std::vector<std::vector<hy_solution_t>> sols_1(ngrid, sols_0);
	std::vector<std::vector<std::vector<hy_solution_t>>> sols(ngrid, sols_1);
	
	std::cout << std::setprecision(15);
	for (int i_r = 0; i_r < ngrid; ++i_r)
	{
		for (int i_t = 0; i_t < ngrid; ++i_t)
		{
			for (int i_u = 0; i_u < ngrid; ++i_u)
			{
				double dr_loc = (i_r-((ngrid-1)/2))*e_r;
				double dt_loc = (i_t-((ngrid-1)/2))*e_t;
				double du_loc = (i_u-((ngrid-1)/2))*e_u;
				std::cout << dr_loc << " (" << dr_loc*HyWall::settings.gasConstant*T_f << "), " << dt_loc << ", " << du_loc << std::endl;
				sols[i_r][i_t][i_u] = solve(du_loc, dt_loc, dr_loc, T_f);
			}
		}
	}
	if (!std::filesystem::exists("deriv")) std::filesystem::create_directory("deriv");
	int nderiv = ngrid-1;
	for (int dr = 0; dr <= nderiv; ++dr)
	{
		for (int dt = 0; dt <= nderiv; ++dt)
		{
			for (int du = 0; du <= nderiv; ++du)
			{
				std::cout << "Compute deriv. (" << du << ", " << dt << ", " << dr << ")" << std::endl;
				auto result = compute_deriv(du, dt, dr, sols, e_u, e_t, e_r);
				
				const double nwmp  = HyWall::settings.rayDim;
				const double ncent = (ngrid-1)/2;
				const auto&  def_data = sols[ncent][ncent][ncent];
				const double ly = std::get<0>(def_data)[nwmp-1];
				const double uf = std::get<1>(def_data)[nwmp-1];
				const double Tf = std::get<2>(def_data)[nwmp-1];
				const double rf = std::get<3>(def_data)[nwmp-1];
				const double mf = std::get<4>(def_data)[nwmp-1];
				const double derivnorm = std::pow(uf, du)*std::pow(Tf, dt)*std::pow(rf, dr);
				const double normalize_M = derivnorm*ly/(mf*uf);
				const double normalize_E = derivnorm*ly*HyWall::settings.fluidPrandtl/(HyWall::settings.fluidCp*mf*Tf);
				
				// const double normalize_M = 1.0;
				// const double normalize_E = 1.0;
				
				const auto& dMdx = std::get<0>(result);
				const auto& dEdx = std::get<1>(result);
				const std::string filename = "deriv/d" + std::to_string(du) + std::to_string(dt) + std::to_string(dr) + "EM.csv";
				std::cout << "output " << filename << std::endl;
				std::ofstream deriv_file(filename);
				for (int i = 0; i < HyWall::settings.rayDim; ++i)
				{
					deriv_file << std::get<0>(sols[0][0][0])[i] << ", " << dMdx[i]*normalize_M << ", " << dEdx[i]*normalize_E << "\n";
				}
			}
		}
	}
	
	HyWall::settings.verboseLevel = tmp;
}

void compute_wm_commutation(hy_probes_t& probes, alt_settings_t& alt_settings, std::vector<double>& q_f)
{
	const std::string filename = alt_settings.commutation_filename;
	std::ifstream myfile(filename);
	std::string line;
	
	std::vector<double> avg_y      (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_mu     (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_mu_t   (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_lam_t  (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_u      (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_T      (HyWall::settings.rayDim, 0.0);
	std::vector<double> avg_rho    (HyWall::settings.rayDim, 0.0);
	
	double dist_f_avg = 0.0;
	double rho_f_avg  = 0.0;
	double mu_f_avg   = 0.0;
	double p_f_avg    = 0.0;
	double u_f_avg    = 0.0;
	double v_f_avg    = 0.0;
	double w_f_avg    = 0.0;
	double T_f_avg    = 0.0;

	int num = 0;

	while (std::getline(myfile, line))
	{
		std::istringstream iss(line);
		double dist, x, y, z, rho, mu, p, u, v, w, T, mu_t, mom_rhs, dp_dx;
		x       = 0.0;
		y       = 0.0;
		z       = 0.0;
		mu_t    = 0.0;
		mom_rhs = 0.0;
		dp_dx   = 0.0;
		iss >> dist;
		iss >> rho;
		iss >> mu;
		iss >> p;
		iss >> u;
		iss >> v;
		iss >> w;
		iss >> T;
		std::cout << line << std::endl;
		
		q_f[0] = p;
		q_f[1] = u;
		q_f[2] = v;
		q_f[3] = w;
		q_f[4] = T;
		q_f[5] = mu_t;
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		
		std::vector<double> v_y      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_mu     (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_mu_t   (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_lam_t  (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_u      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_T      (HyWall::settings.rayDim, 0.0);
		std::vector<double> v_rho    (HyWall::settings.rayDim, 0.0);
		
		HyWall::ProbeSolution(probes.y_prb,    0, &v_y    [0]);
		HyWall::ProbeSolution(probes.u_prb,    0, &v_u    [0]);
		HyWall::ProbeSolution(probes.T_prb,    0, &v_T    [0]);
		HyWall::ProbeSolution(probes.rho_prb,  0, &v_rho  [0]);
		HyWall::ProbeSolution(probes.mu_prb,   0, &v_mu   [0]);
		HyWall::ProbeSolution(probes.mu_t_prb, 0, &v_mu_t [0]);
		HyWall::ProbeSolution(probes.lam_t_prb,0, &v_lam_t[0]);
		
		double alpha = 1.0/(num+1);
		for (int i = 0; i < HyWall::settings.rayDim; ++i)
		{
			avg_y[i]     = (1.0-alpha)*avg_y[i]     + (alpha)*v_y[i];
			avg_mu[i]    = (1.0-alpha)*avg_mu[i]    + (alpha)*v_u[i];
			avg_mu_t[i]  = (1.0-alpha)*avg_mu_t[i]  + (alpha)*v_T[i];
			avg_lam_t[i] = (1.0-alpha)*avg_lam_t[i] + (alpha)*v_rho[i];
			avg_u[i]     = (1.0-alpha)*avg_u[i]     + (alpha)*v_mu[i];
			avg_T[i]     = (1.0-alpha)*avg_T[i]     + (alpha)*v_mu_t[i];
			avg_rho[i]   = (1.0-alpha)*avg_rho[i]   + (alpha)*v_lam_t[i];
		}
		num++;
		
		int nn = HyWall::settings.rayDim;
		
		dist_f_avg = (1.0-alpha)*dist_f_avg + (alpha)*dist;
		rho_f_avg  = (1.0-alpha)*rho_f_avg  + (alpha)*rho;
		mu_f_avg   = (1.0-alpha)*mu_f_avg   + (alpha)*mu;
		p_f_avg    = (1.0-alpha)*p_f_avg    + (alpha)*p;
		u_f_avg    = (1.0-alpha)*u_f_avg    + (alpha)*u;
		v_f_avg    = (1.0-alpha)*v_f_avg    + (alpha)*v;
		w_f_avg    = (1.0-alpha)*w_f_avg    + (alpha)*w;
		T_f_avg    = (1.0-alpha)*T_f_avg    + (alpha)*T;
		
		q_f[0] = p_f_avg;
		q_f[1] = u_f_avg;
		q_f[2] = v_f_avg;
		q_f[3] = w_f_avg;
		q_f[4] = T_f_avg;
		q_f[5] = 0.0;
		HyWall::Solve();
		HyWall::Solve();
		HyWall::Solve();
		
		std::vector<double> msc_y      (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_mu     (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_mu_t   (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_lam_t  (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_u      (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_T      (HyWall::settings.rayDim, 0.0);
		std::vector<double> msc_rho    (HyWall::settings.rayDim, 0.0);
		
		HyWall::ProbeSolution(probes.y_prb,    0, &msc_y[0]);
		HyWall::ProbeSolution(probes.u_prb,    0, &msc_mu[0]);
		HyWall::ProbeSolution(probes.T_prb,    0, &msc_mu_t[0]);
		HyWall::ProbeSolution(probes.rho_prb,  0, &msc_lam_t[0]);
		HyWall::ProbeSolution(probes.mu_prb,   0, &msc_u[0]);
		HyWall::ProbeSolution(probes.mu_t_prb, 0, &msc_T[0]);
		HyWall::ProbeSolution(probes.lam_t_prb,0, &msc_rho[0]);
		
		std::string msc_filename = "msc.csv";
		std::string mfc_filename = "mfc.csv";
		
		std::ofstream m1(msc_filename);
		std::ofstream m2(mfc_filename);
		for (int i = 0; i < nn; ++i)
		{
			m1 << msc_y[i] << ", " << msc_u[i] << ", " << msc_T[i] << ", " << msc_mu[i] << ", " << msc_rho[i] << ", " << msc_mu_t[i] << ", " << msc_lam_t[i] << "\n";
			m2 << avg_y[i] << ", " << avg_u[i] << ", " << avg_T[i] << ", " << avg_mu[i] << ", " << avg_rho[i] << ", " << avg_mu_t[i] << ", " << avg_lam_t[i] << "\n";
		}
	}
}

int main(int argc, char** argv)
{
	std::string filename = "input.ptl";
	if (argc>1) filename = std::string(argv[1]);
	MPI_Init(&argc, &argv);
	bool was_init = false;
	input_data_t input_data;
	HyWall::Initialize(MPI_COMM_WORLD, 4);
	
	alt_settings_t alt_settings;
	std::string out_filename;
	ReadInput(filename, HyWall::settings, input_data, was_init, out_filename, argc, argv, alt_settings);
	
	
	if (std::abs(input_data.rho - input_data.p/(input_data.T*HyWall::settings.gasConstant))>1e-5)
	{
		std::cout << "==============================================" << std::endl;
		std::cout << "WARNING!" << std::endl;
		std::cout << "Failed approximation of ideal gas law:" << std::endl;
		std::cout << "P:   " << input_data.p << std::endl;
		std::cout << "T:   " << input_data.T << std::endl;
		std::cout << "Rho: " << input_data.rho << " (found), " << input_data.p/(input_data.T*HyWall::settings.gasConstant) << " (expected)" << std::endl;
		std::cout << "==============================================" << std::endl;
	}
	if (was_init) return 0;
	HyWall::SetDomainSize(1);
	HyWall::DefineVariables();
	std::vector<double> ff_in(6);
	ff_in[0] = input_data.p;
	ff_in[1] = input_data.u;
	ff_in[2] = input_data.v;
	ff_in[3] = input_data.w;
	ff_in[4] = input_data.T;
	ff_in[5] = input_data.mu_t;
	HyWall::PassFlowfieldVariables(&ff_in[0], 1);
	HyWall::PassVariable("in:distance",    &input_data.dist);
	HyWall::PassVariable("in:x",           &input_data.x);
	HyWall::PassVariable("in:y",           &input_data.y);
	HyWall::PassVariable("in:z",           &input_data.z);
	HyWall::PassVariable("in:rho",         &input_data.rho);
	HyWall::PassVariable("in:mu_lam",      &input_data.mu);
	HyWall::PassVariable("in:momRHS",      &input_data.mom_rhs);
	HyWall::PassVariable("in:dpdx",        &input_data.dp_dx);
	double dummy = 1.0;
	if (HyWall::settings.enableTransitionSensor)
	{
		HyWall::SetTimeStep(1.0);
		HyWall::PassVariable("aux:strain_rate",    &dummy);
		HyWall::PassVariable("aux:sensor_preMult", &dummy);
	}
	
	output_data_t output_data;
	HyWall::PassVariable("out:vorticity",    &output_data.vort);
	HyWall::PassVariable("out:tau",          &output_data.tau);
	HyWall::PassVariable("out:heatflux",     &output_data.qw);
	HyWall::PassVariable("out:failurelevel", &output_data.fail);
	HyWall::Allocate();
	HyWall::Solve();
	
	std::cout << "Finished solve.\n";
	std::cout << "==== INPUTS  ====\n";
	std::cout << " >> dist: " << input_data.dist << "\n";
	std::cout << " >> rho:  " << input_data.rho  << "\n";
	std::cout << " >> mu:   " << input_data.mu   << "\n";
	std::cout << " >> p:    " << input_data.p    << "\n";
	std::cout << " >> u:    " << input_data.u    << "\n";
	std::cout << " >> T:    " << input_data.T    << "\n";
	std::cout << "==== OUTPUTS ====\n";
	std::cout << " >> tau: " << output_data.tau << "\n";
	std::cout << " >> qw:  " << output_data.qw  << "\n";
	
	HyWall::InitProbeIndex();
	hy_probes_t probes;
	HyWall::DefineProbeIndex("sol:d",     &probes.y_prb);
	HyWall::DefineProbeIndex("sol:u",     &probes.u_prb);
	HyWall::DefineProbeIndex("sol:T",     &probes.T_prb);
	HyWall::DefineProbeIndex("sol:mu_t",  &probes.mu_t_prb);
	HyWall::DefineProbeIndex("sol:lam_t", &probes.lam_t_prb);
	HyWall::DefineProbeIndex("sol:rho",   &probes.rho_prb);
	HyWall::DefineProbeIndex("sol:mu",    &probes.mu_prb);
	
	std::vector<double> v_y   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_u   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_T   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_rho (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_mu  (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_mu_t(HyWall::settings.rayDim, 0.0);
	
	HyWall::ProbeSolution(probes.y_prb,    0, &v_y   [0]);
	HyWall::ProbeSolution(probes.u_prb,    0, &v_u   [0]);
	HyWall::ProbeSolution(probes.T_prb,    0, &v_T   [0]);
	HyWall::ProbeSolution(probes.rho_prb,  0, &v_rho [0]);
	HyWall::ProbeSolution(probes.mu_prb,   0, &v_mu  [0]);
	HyWall::ProbeSolution(probes.mu_t_prb, 0, &v_mu_t[0]);
	std::vector<double> prt;
	prt.resize(HyWall::settings.rayDim, 0.0);
	for (int i = 0; i < prt.size(); i++)
	{
		prt[i] = HyCore::GetTurbPrandtl(0, i, HyWall::settings.variablePrandtlT, HyWall::settings.yscaleType);
	}
	std::ofstream out_file(out_filename);
	out_file << std::setprecision(20);
	for (int i = 0; i < v_y.size(); i++)
	{
		out_file << v_y[i] << ", " << v_u[i] << ", " << v_T[i] << ", " << v_rho[i] << ", " << v_mu[i] << ", " << v_mu_t[i] << ", " << prt[i] << "\n";
	}
	
	if (alt_settings.compute_derivatives)
	{
		compute_wm_derivatives(probes, alt_settings, ff_in);
	}
	
	if (alt_settings.compute_commutation)
	{
		compute_wm_commutation(probes, alt_settings, ff_in);
	}
	
	HyWall::Finalize();
	MPI_Finalize();
	return 0;
}
