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

bool streq(const std::string& s1, const std::string& s2) {return s1.compare(s2)==0;}
void ReadInput(const std::string& filename, decltype(HyWall::settings)& settings, input_data_t& in_data, bool& ran_init, std::string& out_file, int argc, char** argv)
{
	PTL::AddUserDefinedFunction("cleandata", cleandata);
	PTL::AddUserDefinedFunction("specialfread", specialfread);
	PTL::AddUserDefinedFunction("suth", suth);
	PTL::PropertyTree input;
	PTL::Interactive i(argc, argv, &input);
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
	
	if (!failed) input.StrictParse();
	     if (streq(mom_eq_str, "allmaras")) {settings.momentumEquationType = HyCore::momentum::allmaras;}
	else if (streq(mom_eq_str, "ODE"))      {settings.momentumEquationType = HyCore::momentum::ODE;}
	else if (!failed) {std::cout << "Invalid value for momentum equation" << std::endl; abort();}

	     if (streq(trb_eq_str, "linear"))    {settings.turbulenceEquationType = HyCore::turbulence::linear;}
	else if (streq(trb_eq_str, "ODE"))       {settings.turbulenceEquationType = HyCore::turbulence::ODE;}
	else if (streq(trb_eq_str, "vanDriest")) {settings.turbulenceEquationType = HyCore::turbulence::vanDriest;}
	else if (!failed) {std::cout << "Invalid value for turbulence equation" << std::endl; abort();}

	     if (streq(eng_eq_str, "croccoBusemann")) {settings.energyEquationType = HyCore::energy::croccoBusemann;}
	else if (streq(eng_eq_str, "ODE"))            {settings.energyEquationType = HyCore::energy::ODE; }
	else if (streq(eng_eq_str, "linear"))         {settings.energyEquationType = HyCore::energy::linear;}
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

int main(int argc, char** argv)
{
	std::string filename = "input.ptl";
	if (argc>1) filename = std::string(argv[1]);
	MPI_Init(&argc, &argv);
	bool was_init = false;
	input_data_t input_data;
	HyWall::Initialize(MPI_COMM_WORLD, 4);
	
	std::string out_filename;
	ReadInput(filename, HyWall::settings, input_data, was_init, out_filename, argc, argv);
	
	
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
	int u_prb, T_prb, mu_t_prb, y_prb, rho_prb, mu_prb;
	HyWall::DefineProbeIndex("sol:d",    &y_prb);
	HyWall::DefineProbeIndex("sol:u",    &u_prb);
	HyWall::DefineProbeIndex("sol:T",    &T_prb);
	HyWall::DefineProbeIndex("sol:mu_t", &mu_t_prb);
	HyWall::DefineProbeIndex("sol:rho",  &rho_prb);
	HyWall::DefineProbeIndex("sol:mu",   &mu_prb);
	
	std::vector<double> v_y   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_u   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_T   (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_rho (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_mu  (HyWall::settings.rayDim, 0.0);
	std::vector<double> v_mu_t(HyWall::settings.rayDim, 0.0);
	
	HyWall::ProbeSolution(y_prb,    0, &v_y   [0]);
	HyWall::ProbeSolution(u_prb,    0, &v_u   [0]);
	HyWall::ProbeSolution(T_prb,    0, &v_T   [0]);
	HyWall::ProbeSolution(rho_prb,  0, &v_rho [0]);
	HyWall::ProbeSolution(mu_prb,   0, &v_mu  [0]);
	HyWall::ProbeSolution(mu_t_prb, 0, &v_mu_t[0]);
	
	std::ofstream out_file(out_filename);
	for (int i = 0; i < v_y.size(); i++)
	{
		out_file << v_y[i] << ", " << v_u[i] << ", " << v_T[i] << ", " << v_rho[i] << ", " << v_mu[i] << ", " << v_mu_t[i] << "\n";
	}
	
	HyWall::Finalize();
	MPI_Finalize();
	return 0;
}
