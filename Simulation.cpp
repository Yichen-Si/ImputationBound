#include "tree2.h"

using namespace std;

class InputParser{
  public:
    InputParser (int &argc, char **argv){
      for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
    }
    std::string getCmdOption(const std::string &option) {
      std::vector<std::string>::iterator itr;
      itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
      if (itr != this->tokens.end() && ++itr != this->tokens.end()){
          return *itr;
      }
      static std::string empty_string("");
      return empty_string;
    }
    bool cmdOptionExists(const std::string &option) const{
      return std::find(this->tokens.begin(), this->tokens.end(), option)
      != this->tokens.end();
    }
  private:
      std::vector <std::string> tokens;
};

int main(int argc, char **argv) {

  int n = 1e3, j = 25, u = 50, h = 0;
  double N0 = 1e6, N1 = 2e3, N2 = 2e4;
  int t0 = 500, t1 = 1500, t2 = 1000;
  double beta = 1.25;
  int seed = 1984;
  int R = 1;
  int config = 0;
  bool tflag = 0;
  double thin = 0.01;

  string ntmp, outfile, prefix;

  InputParser input(argc, argv);

  string stmp = input.getCmdOption("-s");
  if (!stmp.empty())
    seed = stoi(stmp);
  else
    cout << "Use default seed: " << seed << endl;

  ntmp = input.getCmdOption("-thin");
  if (!ntmp.empty())
    thin = stof(ntmp);

  ntmp = input.getCmdOption("-R");
  if (!ntmp.empty())
    R = (int) stof(ntmp);

  ntmp = input.getCmdOption("-n");
  if (!ntmp.empty())
    n = (int) stof(ntmp);

  ntmp = input.getCmdOption("-j");
  if (!ntmp.empty())
    j = stoi(ntmp);

  ntmp = input.getCmdOption("-u");
  if (!ntmp.empty())
    u = stoi(ntmp);

  ntmp = input.getCmdOption("-h");
  if (!ntmp.empty())
    h = stoi(ntmp);

  ntmp = input.getCmdOption("-N0");
  if (!ntmp.empty())
    N0 = stof(ntmp);

  ntmp = input.getCmdOption("-N1");
  if (!ntmp.empty())
    N1 = stof(ntmp);

  ntmp = input.getCmdOption("-N2");
  if (!ntmp.empty())
    N2 = stof(ntmp);

  ntmp = input.getCmdOption("-t0");
  if (!ntmp.empty())
    t0 = (int) stof(ntmp);

  ntmp = input.getCmdOption("-t1");
  if (!ntmp.empty())
    t1 = (int) stof(ntmp);

  ntmp = input.getCmdOption("-t2");
  if (!ntmp.empty())
    t2 = (int) stof(ntmp);

  ntmp = input.getCmdOption("-timeonly");
  if (!ntmp.empty())
    tflag = (bool) stoi(ntmp);

  ntmp = input.getCmdOption("-beta");
  if (!ntmp.empty())
    beta = stod(ntmp);

  ntmp = input.getCmdOption("-count");
  if (!ntmp.empty())
    config = (int) stof(ntmp);

  prefix = input.getCmdOption("-prefix");
  if (prefix.empty())
    prefix = "./test";

  outfile = prefix + "_Prob_n_" + to_string(n) + "_R_" + to_string(R) + "_s_" + to_string(seed) + ".txt";
  cout << "Save to output file: " << outfile << endl;
  cout << "Parameters:" << endl;
  cout << "History type " << h << "\tSample size " << n << endl;
  cout << "Max derived allele count " << j <<  endl;

  clock_t begin = clock();

  srand(seed*1984*1e4+2018);
  std::default_random_engine rng = Rng.getRNG();
  std::bernoulli_distribution bin(thin);

  ofstream ofile;
  ofile.open(outfile);

  for (int r = 0; r < R; ++r) {
    Rng.set_seed(rand());
    // cout << "R " << r << '\n';
    Tree spruce(n, j, u, h, N0, t0, N1, t1, N2, t2, beta);
    // cout << "Build one tree with " << spruce.candynotes.size() << " relative internal branches \n";
    for (auto & v : spruce.candynotes) {
      if (v->size < 10 && !bin(rng)) {continue;}
      double wrong=0.0, nonzero=0.0, overhalf=0.0;
      double one = v->size;
      double tot = v->size;
      double k = v->event;
      double avghat = one, avghatsq = one;
      Node* w = v->parent;
      k = w->event;
      wrong = (k-1.0)/n/(1.0+n);
      while(tot <= one * 2 && k > 2) {
        overhalf += (k-1.0)/n;
        nonzero += (k-1.0)/n;
        avghat +=  (one/tot) * (k-1.0)/n;
        avghatsq += pow(one/tot,2) * (k-1.0)/n;
        tot = w->size;
        w = w->parent;
        k = w->event;
      }
      // cout << "I: " << k << '\n';
      while(k > 2) {
        nonzero += (k-1.0)/n;
        avghat +=  (one/tot) * (k-1.0)/n;
        avghatsq += pow(one/tot,2) * (k-1.0)/n;
        tot = w->size;
        w = w->parent;
        k = w->event;
      }
      // cout << "II: " << k << '\n';
      overhalf /= (1.0+n);
      nonzero /= (1.0+n);
      avghat/= (1.0+n);
      avghatsq/= (1.0+n);

      double r2 = one / (n-one) * pow(1.0 - avghat,2) / (avghatsq - pow(avghat,2.0));
      // cout << n << '\t' << v->size << '\t' << avghat << '\t' << avghatsq << '\t' << r2 << '\t' << v->branch_length_up << '\n';
      ofile << n << '\t' << v->size << '\t' << wrong << '\t' << overhalf << '\t' << nonzero << '\t' << r2 << '\t' << v->branch_length_up << '\n';
    }
  }

  ofile.close();
  return 0;
}


