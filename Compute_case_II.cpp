 #include "Combinatorics.h"

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

  int n = 1e3, maxac, k;
  string ntmp, tfile1, tfile2, outfile, prefix;

  InputParser input(argc, argv);

  ntmp = input.getCmdOption("-n");
  if (!ntmp.empty())
    n = (int) stof(ntmp);

  ntmp = input.getCmdOption("-maxac");
  if (!ntmp.empty())
    maxac = stoi(ntmp);

  tfile1 = input.getCmdOption("-tf1");
  tfile2 = input.getCmdOption("-tf2");

  prefix = input.getCmdOption("-prefix");
  if (prefix.empty())
    prefix = "./test";

  // Read coalescent time
  vector<double> tk1(n, 0);
  vector<double> tk2(n+1, 0);
  ifstream rf;
  rf.open(tfile1, ifstream::in);
  double tk = 0.0;
  int i = 0;
  while (rf >> tk) {
    tk1[i] = tk;
    i++;
  }
  rf.close();
  rf.open(tfile2, ifstream::in);
  i = 0;
  while (rf >> tk) {
    tk2[i] = tk;
    i++;
  }
  if (i < n+1) {
    tk2[i] = tk;
    i++;
  }
  rf.close();
  cout << "Finish reading time\n";

  double **fjkmat;
  fjkmat = FJK(maxac, n);
  cout << "Finishe computing f(j,k;n)\n";

  vector<double> res(maxac+1,0.0);
  vector<double> up(maxac+1,0.0);
  vector<double> down(maxac+1,0.0);

  double err = 1e-6;

  for (k = 2; k <= n - maxac + 1; ++k) {
    double up_in = 0.0, down_in = 0.0;
    for (int d = 2; d < k; ++d) {
      up_in += d * (d-1) * (tk2[d+1] - tk2[k+1]);
      down_in += (d-1) * (tk1[d-1] - tk1[k]);
    }
    for (int j = 1; j <= maxac; ++j) {
      up[j] += up_in * exp(fjkmat[j][k]);
      down[j] += down_in * exp(fjkmat[j][k]);
    }
    if (k > 100 && up[maxac] > 0) {
      double delta_up = up_in * exp(fjkmat[maxac][k]);
      double delta_down = down_in * exp(fjkmat[maxac][k]);
      if (up_in > 0 && delta_up / up[maxac] < err && delta_down / up[maxac] < err) {
        break;
      }
    }
  }

  cout << "Big iteration I, stopped at " << k << '\n';
  cout << "Up: " << up[2] << '\t' << up[maxac] << '\n';
  cout << "Down: " << down[2] << '\t' << down[maxac] << '\n';

  for (k = 2; k <= n - maxac + 1; ++k) {
    double up_in = 0.0;
    for (int k0 = 2; k0 < k; ++k0) {
      double delta_in = 0.0;
      for (int d = k0; d < k; ++d) {
        delta_in += d * (tk2[d+1] - tk2[k+1]);
        // up_in += d * (tk2[d+1] - tk2[k+1]);
      }
      up_in += delta_in;
      if (k0 > 100 && delta_in / up_in < err) {
        break;
      }
    }
    for (int j = 1; j <= maxac; ++j) {
      up[j] += 2.0 * up_in * exp(fjkmat[j][k]);
    }
    double delta_up = 2.0 * up_in * exp(fjkmat[maxac][k]);
    if (k > 100 && delta_up / up[maxac] < err) {
      break;
    }
    if (k % 500 == 0) {
      cout << "Big iteration II, at " << k << "\t Delta_up " << delta_up << "\t Change by " << delta_up / up[maxac] << '\n';
    }
  }

  cout << "Big iteration II, stopped at " << k << '\n';

  res[maxac] = up[maxac] / down[maxac] / ((double) (n*(n+1)));
  cout << maxac << '\t' << up[maxac] << '\t' << down[maxac] << '\t' << res[maxac]  << '\n';

  for (int ac = (maxac-1); ac > 0; --ac) {
    k = n - ac + 1;
    double up_in = 0.0, down_in = 0.0;
    for (int d = 2; d < k; ++d) {
      up_in += d * (d-1) * (tk2[d+1] - tk2[k+1]);
      down_in += (d-1) * (tk1[d-1] - tk1[k]);
    }
    for (int j = 1; j <= ac; ++j) {
      up[j] += up_in * exp(fjkmat[j][k]);
      down[j] += down_in * exp(fjkmat[j][k]);
    }

    up_in = 0.0;
    for (int k0 = 2; k0 < k; ++k0) {
      for (int d = k0; d < k; ++d) {
        up_in += d * (tk2[d+1] - tk2[k+1]);
      }
    }
    for (int j = 1; j <= ac; ++j) {
      up[j] += 2.0 * up_in * exp(fjkmat[j][k]);
    }

    res[ac] = up[ac] / down[ac] / ((double) (n*(n+1)));
  cout << ac << '\t' << up[ac] << '\t' << down[ac] << '\t' << res[ac]  << '\n';
  }

  outfile = prefix + "_ProbPartial.txt";
  ofstream sfile;
  sfile.open(outfile);
  for (int j = 0; j <= maxac; ++j) {
    sfile << j << '\t' << res[j] << endl;
  }
  sfile.close();

}
