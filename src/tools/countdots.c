#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdbool.h>
#include <glib.h>

typedef struct
{
    char *node1;
    char *node2;
} Edge;

typedef struct
{
    GHashTable *node_set;
    GHashTable *edge_set;
    int num_nodes;
    long num_edges;
} Graph;

void add_node(Graph *graph, const char *node)
{
    if (!g_hash_table_contains(graph->node_set, node))
    {
        g_hash_table_add(graph->node_set, g_strdup(node));
        graph->num_nodes++;
    }
}

void add_edge(Graph *graph, const char *node1, const char *node2)
{
    // Sort the nodes before adding the edge
    const char *sorted_node1;
    const char *sorted_node2;
    if (strcmp(node1, node2) < 0)
    {
        sorted_node1 = node1;
        sorted_node2 = node2;
    }
    else
    {
        sorted_node1 = node2;
        sorted_node2 = node1;
    }

    // Create the edge key dynamically
    char *edge_key = g_strdup_printf("%s-%s", sorted_node1, sorted_node2);

    // Check if the edge already exists in the set
    if (g_hash_table_contains(graph->edge_set, edge_key)) {
        g_free(edge_key);
        return;
    }
    // Add the edge to the set
    g_hash_table_add(graph->edge_set, edge_key);
    graph->num_edges++;
}

char *extract_node_name(const char *input, bool is_second_node)
{
    if (!input) return NULL;

    // Trim leading whitespace
    while (isspace((unsigned char)*input))
    {
        input++;
    }

    if (*input == '\0') return NULL;

    char *result = NULL;
    if (*input == '\"')
    {
        // Enclosed in double quotes
        input++; // Skip opening quote
        const char *end = strchr(input, '\"');
        if (end)
        {
            size_t len = end - input;
            result = (char *)malloc((len + 1) * sizeof(char));
            if (result)
            {
                strncpy(result, input, len);
                result[len] = '\0';
            }
        }
        else
        {
            // Unclosed quote, extract everything
            result = strdup(input);
        }
    }
    else
    {
        // Not enclosed in double quotes
        if (is_second_node)
        {
            // For the second node, read until whitespace, ';', '[', or '\0'
            const char *p = input;
            while (*p && !isspace((unsigned char)*p) && *p != ';' && *p != '[')
            {
                p++;
            }
            size_t len = p - input;
            result = (char *)malloc((len + 1) * sizeof(char));
            if (result)
            {
                strncpy(result, input, len);
                result[len] = '\0';
            }
        }
        else
        {
            // For the first node, it's the entire remaining string (trimmed)
            const char *end = input + strlen(input) - 1;
            while (end > input && isspace((unsigned char)*end))
            {
                end--;
            }
            size_t len = end - input + 1;
            result = (char *)malloc((len + 1) * sizeof(char));
            if (result)
            {
                strncpy(result, input, len);
                result[len] = '\0';
            }
        }
    }
    return result;
}

void count_nodes_edges(const char *dot_file, Graph *graph)
{
    FILE *file = fopen(dot_file, "r");
    if (file == NULL)
    {
        printf("Failed to open file: %s\n", dot_file);
        exit(1);
    }

    // Create the sets with key destroy notification to prevent memory leaks
    graph->edge_set = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    graph->node_set = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

    char line[65536];
    while (fgets(line, sizeof(line), file))
    {
        char *node1 = NULL;
        char *node2 = NULL;

        // Check if the line contains an edge connection
        char *sep = strstr(line, " -- ");
        if (sep)
        {
            *sep = '\0';
            node1 = extract_node_name(line, false);
            node2 = extract_node_name(sep + 4, true);
        }

        if (node1 != NULL && node2 != NULL)
        {
            add_edge(graph, node1, node2);
            add_node(graph, node1);
            add_node(graph, node2);
        }

        if (node1 != NULL) free(node1);
        if (node2 != NULL) free(node2);
    }

    fclose(file);
}

void free_graph(Graph *graph)
{
    if (graph)
    {
        if (graph->edge_set)
        {
            g_hash_table_destroy(graph->edge_set);
        }
        if (graph->node_set)
        {
            g_hash_table_destroy(graph->node_set);
        }
        free(graph);
    }
}

Graph *process_file(const char *dot_file, int print_labels)
{
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    if (graph == NULL)
    {
        printf("Memory allocation failed\n");
        return NULL;
    }
    graph->num_nodes = 0;
    graph->num_edges = 0;
    graph->edge_set = NULL;
    graph->node_set = NULL;

    count_nodes_edges(dot_file, graph);

    if (print_labels)
    {
        printf("File\tNodes\tEdges\n");
    }

    char *filename = strrchr(dot_file, '/');
    if (filename == NULL)
    {
        filename = (char *)dot_file;
    }
    else
    {
        filename++;
    }

    printf("%s\t%d\t%ld\n", filename, graph->num_nodes, graph->num_edges);
    return graph;
}

void process_folder(const char *folder)
{
    int print_labels = 1;
    int total_nodes = 0;
    int total_edges = 0;

    char command[256];
    sprintf(command, "find %s -name '*.dot'", folder);
    FILE *pipe = popen(command, "r");
    if (pipe == NULL)
    {
        printf("Failed to execute command: %s\n", command);
        exit(1);
    }
    char dot_file[256];
    while (fgets(dot_file, sizeof(dot_file), pipe))
    {
        dot_file[strcspn(dot_file, "\n")] = 0;

        Graph *graph = process_file(dot_file, print_labels);
        print_labels = 0;

        if (graph != NULL)
        {
            total_nodes += graph->num_nodes;
            total_edges += graph->num_edges;
            free_graph(graph);
        }
    }

    pclose(pipe);

    if (total_nodes > 0 || total_edges > 0)
    {
        printf("\nSummary - Total nodes: %d, Total edges: %d\n", total_nodes, total_edges);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: ./countdots <dot_file_or_folder>\n");
        return 1;
    }

    char *path = argv[1];

    if (access(path, F_OK) != 0)
    {
        printf("Invalid file or folder path.\n");
        return 1;
    }

    if (access(path, R_OK) != 0)
    {
        printf("Permission denied: %s\n", path);
        return 1;
    }
    struct stat path_stat;
    if (stat(path, &path_stat) == 0)
    {
        if (S_ISREG(path_stat.st_mode))
        {
            Graph *graph = process_file(path, 1);
            if (graph != NULL)
            {
                free_graph(graph);
            }
        }
        else if (S_ISDIR(path_stat.st_mode))
        {
            process_folder(path);
        }
        else
        {
            printf("Invalid file or folder path.\n");
            return 1;
        }
    }
    else
    {
        printf("Failed to get file or directory status.\n");
        return 1;
    }
    return 0;
}
