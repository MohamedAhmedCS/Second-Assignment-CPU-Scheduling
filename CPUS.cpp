// Mohamed Ahmed
// CSC 139
// Assignment 3
// File: CPUS.cpp
// CPU Scheduling Algorithms
// Section: 03
// Date: 4/28/2024
// Profesor: Fernando Cantillo

// Description: This program reads a file containing process data and a scheduling algorithm, then schedules 
// the processes using the one of the algorithms (Round Robin, Shortest Job First, Priority Schedulin without 
// Preemption, Priority Scheduling with Preemption). The scheduling results are written to an output file.

// The program is run from the command line using : g++ CPUS.cpp -o a.out
// with the input file name and output file name as arguments: ./a.out input.txt output.txt

// The input file should have the following format:
// 1. The first line contains the name of the scheduling algorithm (RR, SJF, PR_noPREMP, PR_PREMP).
// 2. The second line contains the quantum time for Round Robin (RR) algorithm.
// 3. The third line contains the number of processes.
// 4. Each subsequent line contains the process data in the following format:
//    <PID> <Arrival Time> <Burst Time> <Priority>
//    where PID is the process ID, Arrival Time is the time at which the process arrives, Burst Time 
//    is the CPU burst time, and Priority is the priority of the process.

// The output file will contain the scheduling results in the following format:
// 1. The first line contains the name of the scheduling algorithm.
// 2. Each subsequent line contains the scheduled process data in the following format:
//    <Start Time> <PID>
//    where Start Time is the time at which the process starts executing and PID is the process ID.
// 3. The last line contains the average waiting time for all processes.

// The program will display an error message if the input file cannot be opened, the number of processes is invalid,
// or the process data is invalid. The program will also display an error message if the output file cannot be opened.
// The program will display an error message if the scheduling algorithm is not supported.

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <string>
#include <stdexcept>
#include <cstdlib>

using namespace std;

// Class representing a process
class Process {
public:
    int pid;         // Process ID
    int arrT;        // Arrival time
    int burstT;      // CPU burst time
    int priority;    // Priority (lower number indicates higher priority)
    int remainingT;  // Remaining burst time for preemptive algorithms
    int waitT;       // Waiting time
};

// Function to read input data from a file
void readInput(const string& filename, vector<Process>& processes, string& algorithm, int& quantum) {
    ifstream inFile(filename);
    if (!inFile.is_open()) {
        throw runtime_error("Error opening input file: " + filename);
    }

    string line;
    getline(inFile, line);  // First line contains the algorithm name
    stringstream ss(line);
    ss >> algorithm;

    // Extract quantum if Round Robin (RR) is selected
    if (algorithm == "RR") {
        ss >> quantum;
    }

    // Read the number of processes
    int numProcesses;
    if (!(inFile >> numProcesses) || numProcesses <= 0) {
        throw runtime_error("Invalid number of processes specified");
    }

    processes.resize(numProcesses);
    for (int i = 0; i < numProcesses; ++i) {
        if (!(inFile >> processes[i].pid >> processes[i].arrT >> processes[i].burstT >> processes[i].priority) ||
            processes[i].pid < 0 || processes[i].arrT < 0 || processes[i].burstT <= 0 || processes[i].priority < 0) {
            throw runtime_error("Invalid process data at line " + to_string(i + 2));
        }
        processes[i].remainingT = processes[i].burstT;
        processes[i].waitT = 0;
    }
    inFile.close();
}

// Function to write the output scheduling results to a file
void writeOutput(const string& filename, const string& algorithm, const vector<pair<int, int>>& schedule, double avgWaitTime) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error opening output file: " << filename << endl;
        exit(1);
    }
    outFile << algorithm << "\n";
    for (const auto& entry : schedule) {
        outFile << entry.first << "   " << entry.second << "\n";
    }
    outFile << "AVG Waiting Time: " << fixed << setprecision(2) << avgWaitTime << "\n";
    outFile.close();
}

// Round Robin (RR) Scheduling
void RoundRobin(const vector<Process>& processes, int quantum, vector<pair<int, int>>& schedule, double& avgWaitTime) {
    queue<int> q;
    vector<int> remainingTime(processes.size());
    vector<int> waitTimes(processes.size(), 0);
    int currentTime = 0;

    // Initialize remaining burst times and check for initial arrivals
    for (int i = 0; i < processes.size(); i++) {
        remainingTime[i] = processes[i].burstT;
        if (processes[i].arrT <= currentTime) {
            q.push(i);
        }
    }

    // Main scheduling loop for Round Robin
    while (!q.empty()) {
        int idx = q.front();  // Get the next process in the queue
        q.pop();

        // Determine the time slice to use for this process
        int startTime = currentTime;
        int timeSlice = min(quantum, remainingTime[idx]);
        currentTime += timeSlice;
        remainingTime[idx] -= timeSlice;

        // Record the schedule
        schedule.push_back({startTime, processes[idx].pid});

        // Add new arrivals to the queue
        for (int i = 0; i < processes.size(); i++) {
            if (processes[i].arrT > startTime && processes[i].arrT <= currentTime && remainingTime[i] == processes[i].burstT) {
                q.push(i);
            }
        }

        // Requeue the current process if it's not done
        if (remainingTime[idx] > 0) {
            q.push(idx);
        } else {
            waitTimes[idx] = currentTime - processes[idx].burstT - processes[idx].arrT;
        }
    }

    // Calculate the average waiting time
    avgWaitTime = accumulate(waitTimes.begin(), waitTimes.end(), 0.0) / processes.size();
}

// Shortest Job First (SJF) Scheduling
void ShortestJobFirst(const vector<Process>& processes, vector<pair<int, int>>& schedule, double& avgWaitTime) {
    vector<int> waitTimes(processes.size(), 0);
    vector<Process> procs = processes;

    // Create a priority queue to prioritize jobs with shorter burst times
    auto cmp = [&procs](int a, int b) {
        return procs[a].burstT > procs[b].burstT || (procs[a].burstT == procs[b].burstT && procs[a].pid > procs[b].pid);
    };
    priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
    int currentTime = 0;
    size_t index = 0;

    // Add arriving processes to the priority queue or adjust the current time
    while (index < procs.size() || !pq.empty()) {
        while (index < procs.size() && procs[index].arrT <= currentTime) {
            pq.push(index++);
        }

        if (pq.empty()) {
            currentTime = procs[index].arrT;
            continue;
        }

        int idx = pq.top();
        pq.pop();
        schedule.push_back({currentTime, procs[idx].pid});
        currentTime += procs[idx].burstT;
        waitTimes[idx] = currentTime - procs[idx].burstT - procs[idx].arrT;
    }

    // Calculate the average waiting time
    avgWaitTime = accumulate(waitTimes.begin(), waitTimes.end(), 0.0) / waitTimes.size();
}

// Priority Scheduling without Preemption
void PriorityNoPreempt(const vector<Process>& processes, vector<pair<int, int>>& schedule, double& avgWaitTime) {
    vector<int> waitTimes(processes.size(), 0);
    vector<Process> procs = processes;

    // Sort processes by arrival time and priority
    sort(procs.begin(), procs.end(), [](const Process& a, const Process& b) {
        if (a.arrT != b.arrT) return a.arrT < b.arrT;
        return a.priority < b.priority;
    });

    // Create a priority queue based on process priority
    auto cmp = [&procs](int a, int b) {
        return procs[a].priority > procs[b].priority || (procs[a].priority == procs[b].priority && procs[a].pid > procs[b].pid);
    };
    priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
    int currentTime = 0;
    size_t index = 0;

    // Process all initially available processes
    while (index < procs.size() && procs[index].arrT <= currentTime) {
        pq.push(index++);
    }

    // Main scheduling loop
    while (!pq.empty()) {
        int idx = pq.top();
        pq.pop();

        // Advance time if the process has not yet arrived
        if (currentTime < procs[idx].arrT) {
            currentTime = procs[idx].arrT;
        }

        // Schedule the selected process
        schedule.push_back({currentTime, procs[idx].pid});
        waitTimes[idx] = currentTime - procs[idx].arrT;
        currentTime += procs[idx].burstT;

        // Add new processes arriving by the current time to the queue
        while (index < procs.size() && procs[index].arrT <= currentTime) {
            pq.push(index++);
        }
    }

    // Calculate the average waiting time
    avgWaitTime = accumulate(waitTimes.begin(), waitTimes.end(), 0.0) / processes.size();
}

// Main function: Reads input data, runs the scheduling algorithms, and writes the results to a file
int main(int argc, char* argv[]) {
    if (argc != 3) {
cerr << "To use this program, please include your input file name and the output file name "
     << "like this: " << argv[0] << " <input file> <output file>." << endl;
        return 1;
    }

    vector<Process> processes;
    string algorithm;
    int quantum = 0;
    string inputFilename = argv[1];
    string outputFilename = argv[2];

    try {
        readInput(inputFilename, processes, algorithm, quantum);
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return 1;
    }

    vector<pair<int, int>> schedule;
    double avgWaitTime = 0.0;

    // Choose the appropriate scheduling algorithm
    if (algorithm == "RR") {
        RoundRobin(processes, quantum, schedule, avgWaitTime);
    } else if (algorithm == "SJF") {
        ShortestJobFirst(processes, schedule, avgWaitTime);
    } else if (algorithm == "PR_noPREMP") {
        PriorityNoPreempt(processes, schedule, avgWaitTime);
    } else {
        cerr << "Unsupported scheduling algorithm: " << algorithm << endl;
        return 1;
    }

    writeOutput(outputFilename, algorithm, schedule, avgWaitTime);
    return 0;
}
