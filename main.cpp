#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <bitset>
#include <random>
#include <chrono>
#include <fstream>
#include <utility>

#define POPULATION 20

std::ofstream fout("Evolutie.txt");

class CCHROMOSOME{

private:

    static int m_nextID;

    int m_ID;
    double m_realValue;
    std::string m_binaryRepresentation;
    double m_fitnessValue;

public:

    CCHROMOSOME(int id, int p_lowerBound, int p_upperBound, int p_precision, int a, int b, int c){
        m_ID = id;
        m_realValue = CCHROMOSOME::random_real(p_lowerBound, p_upperBound, p_precision);
        m_binaryRepresentation = CCHROMOSOME::realToBinary(m_realValue, p_lowerBound,
                                                           p_upperBound, p_precision);
        m_fitnessValue = CCHROMOSOME::function_value(a, b, c, m_realValue);
    }


    [[nodiscard]] int getID() const {
        return m_ID;
    }

    [[nodiscard]] double getRealVal() const{
        return m_realValue;
    }

    [[nodiscard]] std::string getBinRep() const{

        return m_binaryRepresentation;
    }

    [[nodiscard]] double getFitness() const{

        return m_fitnessValue;
    }

    void setID(int new_val){
        this->m_ID = new_val;
    }

    void setBinRep(std::string new_val){
        this->m_binaryRepresentation = std::move(new_val);
    }

    void setRealVal(double new_val){
        this->m_realValue = new_val;
    }

    void setFitnessVal(double new_val){
        this->m_fitnessValue = new_val;
    }

private:

    std::string realToBinary(double realNumber, double lowerBound, double upperBound, int precision) {
        int length = 22;
        int numberOfIntervals = pow(2, length);
        double intervalSize = (upperBound - lowerBound) / numberOfIntervals;
        int intervalIndex = (realNumber - lowerBound) / intervalSize;

        std::string binaryRepresentation = std::bitset<22>(intervalIndex).to_string();
        return binaryRepresentation;
    }

    double random_real(double lower, double upper, int precision) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        static std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with current time
        std::uniform_real_distribution<double> dis(lower, upper);
        double result = dis(gen);
        double factor = pow(10, precision);
        return (result * factor) / factor;
    }

    double function_value(int a, int b, int c, double x)
    {
        return (a*(x*x) + b * x + c);
    }

};

int CCHROMOSOME::m_nextID = 1;

int population = 20;
int lower_bound = -1;
int upper_bound = 2;
int precision = 6;
int a = -1;
int b = 1;
int c = 2;
float length = std::ceil(std::log2((upper_bound - lower_bound) * std::pow(10, precision)));
float discretization = (upper_bound - lower_bound) / (std::pow(2, length));
double F;

double crossoverProb = 0.25;
double mutationProb = 0.01;

void storeAndPrintInitialData(int population, std::vector<CCHROMOSOME>& vec, bool print){
    if(print)
        fout<<"Populatia initiala\n";

    for(int i =0 ; i < population; i++)
    {

        CCHROMOSOME C = CCHROMOSOME(i+1, lower_bound, upper_bound, precision, a, b, c);
        vec.emplace_back(C);
        F += C.getFitness();
        if(print)
            fout<<" "<<C.getID()<<": "<< C.getBinRep()
                     <<" x = "<<  std::fixed << std::setprecision(precision)<< C.getRealVal()
                     <<" f = "<< C.getFitness()<<"\n";

    }
    if(print)
        fout<<"|-----------------------------------------------------|\n";
    fout<<std::flush;

}

void selectionProbability(std::vector<CCHROMOSOME>& vec, bool print){
    if(print)
        fout<<"Probabilitati selectie\n";
    for(const auto& ch : vec)
    {
        double probability = ch.getFitness() / F;
        if(print)
            fout<<" cromozom "<<ch.getID()<<" probability "<< probability <<"\n";
    }
    if(print)
        fout<<"|--------------------------------|\n";
    fout<<std::flush;
}


double random_real(double lower, double upper) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<double> dis(lower, upper);
    double result = dis(gen);
    return result ;
}

double function_value(double x)
{
    return (a*(x*x) + b * x + c);
}

int binSearch(const std::vector<double>& intervals, double num) {
    int low = 0;
    int high = intervals.size() - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;

        // Check if the number falls within the current interval
        if (intervals[mid] <= num && num < intervals[mid + 1]) {
            return mid;
        }
            // If the number is less than the current interval, search in the left half
        else if (num < intervals[mid]) {
            high = mid - 1;
        }
            // If the number is greater than the current interval, search in the right half
        else {
            low = mid + 1;
        }
    }

    // If the number is not within any interval, return -1 or handle as appropriate
    return -1;
}

std::vector<int> probabilityIntervals(std::vector<CCHROMOSOME>& vChromosome, int population, bool print){
    if(print)
        fout<<"Intervale probabilitati selectie\n";

    double sum = 0;
    std::vector<std::pair<double, double>> intervals;
    std::vector<double> values;
    for(int i = 0 ; i < vChromosome.size() + 1; i++)
    {
        double elem = sum / F;
        if(print)
            fout<<elem<<" ";

        values.emplace_back(elem);
        sum += vChromosome[i].getFitness();

    }
    if(print)
        fout<<"\n|--------------------------------|\n";
    fout<<std::flush;

    for(int i = 0; i < values.size(); i++){
        if(i == values.size() - 1)
            break;
        else{
            intervals.emplace_back(std::make_pair(values[i], values[i+1]));
        }
    }

    std::vector<int> selected_chromosoms;
    for(int i = 0 ; i < population; i++)
    {
        double u = random_real(0, 1);
        int u_chomosom = binSearch(values, u);
        selected_chromosoms.emplace_back(u_chomosom+1);
        if(print)
         fout<<"u = " << u<<"  selectam cromozomul " << u_chomosom + 1<<"\n";

    }

    return selected_chromosoms;
}


int randomIndex(int maxIndex) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 gen(seed);
    std::uniform_int_distribution<int> dis(0, maxIndex - 1);
    return dis(gen);
}

std::pair<std::string, std::string> CrossOver(const std::string& str1, const std::string& str2, int index) {
    // Check if index is within the bounds of the strings
    if (index < 0 || index >= static_cast<int>(str1.length()) || index >= static_cast<int>(str2.length())) {
        // Return the input strings unchanged if index is out of bounds
        return std::make_pair(str1, str2);
    }

    // Split the strings at the specified index
    std::string part1_str1 = str1.substr(0, index);
    std::string part2_str1 = str1.substr(index);
    std::string part1_str2 = str2.substr(0, index);
    std::string part2_str2 = str2.substr(index);

    // Concatenate the parts in the desired order
    std::string new_str1 = part1_str1 + part2_str2;
    std::string new_str2 = part1_str2 + part2_str1;

    return std::make_pair(new_str1, new_str2);
}

double binaryToReal(const std::string& binaryString) {
    int length = binaryString.length();
    int numberOfIntervals = pow(2, length);
    double intervalSize = ((double)upper_bound - (double)lower_bound) / numberOfIntervals;

    int intervalIndex = std::stoi(binaryString, nullptr, 2);

    double realNumber = lower_bound + intervalIndex * intervalSize;

    return realNumber;
}

std::string flipBits(std::string& binaryNumber, std::mt19937& gen) {
    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::cout<<dis(gen)<<"\n";
    for (char& bit : binaryNumber) {
        if (dis(gen) < mutationProb) {
            bit = (bit == '0') ? '1' : '0'; // Flip the bit
        }
    }
    return binaryNumber;
}

int main1() {

    std::vector<CCHROMOSOME> vChromosome;
    storeAndPrintInitialData(population, vChromosome, true);
    selectionProbability(vChromosome, true);
    std::vector<int> selected_chromosoms = probabilityIntervals(vChromosome, population, true);

    std::vector<CCHROMOSOME> newChromosome;
    fout<<"Dupa selectie\n";
    for(int i = 0; i < population; i++){
        for(auto& C : vChromosome)
            if(C.getID() == selected_chromosoms[i]){
                newChromosome.emplace_back(C);
                fout<<" "<<i+1<<": "<< C.getBinRep()
                         <<" x = "<<  std::fixed << std::setprecision(precision)<< C.getRealVal()
                         <<" f = "<< C.getFitness()<<" ID: "<<C.getID()<<"\n";

            }
    }
    fout<<"|-----------------------------------------------------|\n";
    fout<<std::flush;

    std::cout<<"NewChromosome\n";
    int j = 1;
    int id;
    for(auto& C : newChromosome){
        id = C.getID();
        C.setID(j);
            std::cout<<" "<<C.getID()<<": "<< C.getBinRep()
                <<" x = "<<  std::fixed << std::setprecision(precision)<< C.getRealVal()
                <<" f = "<< C.getFitness()<<" Old ID: "<<id<<"\n";
                j++;
        }
    std::cout<<"|-----------------------------------------------------|\n";
    std::cout<<std::flush;

    fout<<"Probabilitatea de incrucisare:"<< crossoverProb<<"\n";
    int i = 1;
    std::vector<CCHROMOSOME> selectedForCrossover;
    //std::vector<std::pair<int, CCHROMOSOME>> SelectedForCrossover;
    for(const auto& C : newChromosome)
    {
        double u = random_real(0, 1);
        fout<<i<<": "<<C.getBinRep() <<" u = "<<u;
        if(u < crossoverProb) {
            selectedForCrossover.emplace_back(C);
            //SelectedForCrossover.emplace_back(std::make_pair(i, C));
            fout << " < " << crossoverProb << " participa\n";
        }
        else
            fout<<"\n";
        i++;
    }

    std::cout<<"selectedForCrossover: \n";
    int sid = 1;
    for(const auto& C : selectedForCrossover) {

        std::cout << " " << C.getID() << ": " << C.getBinRep()
                  << " x = " << std::fixed << std::setprecision(precision) << C.getRealVal()
                  << " f = " << C.getFitness() << " Old ID: " << sid << "\n";
        sid++;
    }

    ///--IN FISIER--//

    if (selectedForCrossover.size() % 2 == 0) {
        for (int k = 0; k < selectedForCrossover.size(); k += 2) {
            int crossoverPoint = randomIndex(population);
            std::string bin1 = selectedForCrossover[k].getBinRep();
            std::string bin2 = selectedForCrossover[k + 1].getBinRep();

            std::pair<std::string, std::string> swappedStrings = CrossOver(bin1, bin2, crossoverPoint);
            fout << "Recombinarea dintre cromozomul "
                      << selectedForCrossover[k].getID() << " cu cromozomul "
                      << selectedForCrossover[k + 1].getID() << ":\n"
                      << bin1<<" "
                      << bin2
                      << " punct " << crossoverPoint<<"\n"
                      << "Rezultat " << swappedStrings.first <<" "
                                     << swappedStrings.second<<"\n";

            double x1 = binaryToReal(swappedStrings.first);
            selectedForCrossover[k].setBinRep(swappedStrings.first);
            selectedForCrossover[k].setRealVal(x1);
            selectedForCrossover[k].setFitnessVal(function_value(x1));

            double x2 = binaryToReal(swappedStrings.second);
            selectedForCrossover[k + 1].setBinRep(swappedStrings.second);
            selectedForCrossover[k + 1].setRealVal(x2);
            selectedForCrossover[k].setFitnessVal(function_value(x2));

        }
    } else {
        for (int k = 0; k < selectedForCrossover.size() - 1; k += 2) {
            int crossoverPoint = randomIndex(population);

            std::string bin1 = selectedForCrossover[k].getBinRep();
            std::string bin2 = selectedForCrossover[k + 1].getBinRep();

            std::pair<std::string, std::string> swappedStrings = CrossOver(bin1, bin2, crossoverPoint);
            fout << "Recombinarea dintre cromozomul "
                 << selectedForCrossover[k].getID() << " cu cromozomul "
                 << selectedForCrossover[k + 1].getID() << ":\n"
                 << selectedForCrossover[k].getBinRep()<<" "
                 << selectedForCrossover[k + 1].getBinRep()
                 << " punct " << crossoverPoint<<"\n"
                    << "Rezultat " << swappedStrings.first <<" "
                    << swappedStrings.second<<"\n";

            double x1 = binaryToReal(swappedStrings.first);
            selectedForCrossover[k].setBinRep(swappedStrings.first);
            selectedForCrossover[k].setRealVal(x1);
            selectedForCrossover[k].setFitnessVal(function_value(x1));

            double x2 = binaryToReal(swappedStrings.second);
            selectedForCrossover[k + 1].setBinRep(swappedStrings.second);
            selectedForCrossover[k + 1].setRealVal(x2);
            selectedForCrossover[k].setFitnessVal(function_value(x2));
        }
    }

    for(const auto& SC  : selectedForCrossover){
        for(auto& C : newChromosome){
            if(SC.getID() == C.getID()){
                C.setBinRep(SC.getBinRep());
                C.setRealVal(SC.getRealVal());
                C.setFitnessVal(SC.getFitness());
            }
        }
    }

    fout<<"Dupa recombinare: \n";
    for(auto& C : newChromosome){
        fout << " " << C.getID() << ": " << C.getBinRep()
                  << " x = " << std::fixed << std::setprecision(precision) << C.getRealVal()
                  << " f = " << C.getFitness() << "\n";
    }

///-------------------------MUTATIE-----------------------------///

    fout<<"Probabilitate de mutatie pentru fiecare gena " << mutationProb
        << "\nAu fost modificati cromozomii:\n";
    double MutationProb = 0.01;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    for(auto& C : newChromosome){
        std::string bin = C.getBinRep();
        //srand(static_cast<unsigned int>(time(nullptr)));
        std::string flipedBin = flipBits(bin, gen);
        if(C.getBinRep() != flipedBin)
        {
            fout<<C.getID()<<"\n";
            C.setBinRep(flipedBin);
            double x = binaryToReal(flipedBin);
            C.setRealVal(x);
            C.setFitnessVal(x);
        }

    }
    fout<<"Dupa mutatie:\n";
    for(auto& C : newChromosome) {
        fout << " " << C.getID() << ": " << C.getBinRep()
             << " x = " << std::fixed << std::setprecision(precision) << C.getRealVal()
             << " f = " << C.getFitness() << "\n";
    }


    return 0;
}


int main() {

    const int iterations = 50;
    for(int i = 0; i < iterations; i++)
    {
        if(i == 0)
            main1();
        else {
            fout<<i<<"\n";
        }
    }
    return 0;
}

